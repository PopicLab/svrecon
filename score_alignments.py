import argparse
import mappy
import pysam
import os
import tempfile

from collections import defaultdict, Counter
from typing import Dict, List

from insilicosv.simulate import SVSimulator


class AlignScorer(object):
    def __init__(self, ref_sequence, sample_sequence, callset_vcf, output_dir, buffer):
        """

        :param ref_sequence: Sequence file for reference genome
        :param sample_sequence: Sequence file for sample genome
        :param callset_vcf: File containing SV callset
        :param output_dir: Output directory for auxiliary files
        :param buffer: Context buffer length
        """
        self.variants = None
        self.vcf_header = None
        self.read_vcf(callset_vcf)
        self.aligner = mappy.Aligner(fn_idx_in=sample_sequence, preset='map-pb')
        print(f"Loaded aligner with sequences {self.aligner.seq_names}")
        self.ref_fasta = pysam.FastaFile(ref_sequence)
        self.sample_sequence = sample_sequence
        self.working_dir = output_dir
        self.buffer = buffer

    def read_vcf(self, vcf_path: str):
        """
        Read callset .vcf and group SVs by SVID info field
        :param vcf_path: path to .vcf file
        """
        grouped_variants = defaultdict(list)
        vcf_in = pysam.VariantFile(vcf_path)

        # Iterate over all records in the VCF file
        for rec in vcf_in.fetch():
            # Check the INFO field for 'SVID'
            svid = rec.info.get('SVID')
            if svid:
                grouped_variants[svid].append(rec)

        self.variants = grouped_variants
        self.vcf_header = vcf_in.header

    def simulate_subsequences(self, svid: str) -> List[Dict]:
        """
        Execute operations defined by called SV on the reference genome
        :param svid: SVID key for SV
        :return: List of sequence dictionaries, each element containing "chrom", "sequence", and "location"

        For SVs without dispersions, the list should only contain a single element. For dispersions,
        the list should include the evidence of the SV at the source and the target.
        """
        records = self.variants[svid]

        sv_type = records[0].info['SVTYPE']

        chrom = records[0].chrom

        start = min([rec.start for rec in records])
        stop = max([rec.stop for rec in records])

        start_buffer = min(self.buffer, start)
        stop_buffer = min(stop + self.buffer, self.ref_fasta.get_reference_length(chrom))

        orig_sequence = self.ref_fasta.fetch(reference=chrom,
                                             start=start - start_buffer,
                                             end=stop + stop_buffer)

        new_sequence = list(orig_sequence)

        subsequences = []  # multiple subsequences may occur if there is a dispersion

        # Create the subsequence corresponding to the SV
        if sv_type in ['INVdel']:
            # No dispersions. Directly apply insilicoSV operations.
            # We're implementing this with O(N) operations for simplicity

            delete_placeholder = ''

            for rec in records:
                offset_start = rec.start - start - start_buffer
                offset_stop = rec.stop - start - start_buffer
                if rec.info['OP_TYPE'] == 'CUT':
                    new_sequence[offset_start:offset_stop] = [delete_placeholder] * (offset_stop - offset_start)
                if rec.info['OP_TYPE'] == 'INV':
                    new_sequence[offset_start:offset_stop] = new_sequence[offset_start:offset_stop:-1]

            new_sequence = ''.join(new_sequence)

            subsequence = {
                'sequence': new_sequence,
                'chrom': chrom,
                'location': start
            }

            subsequences.append(subsequence)

        return subsequences

    def match_subsequence(self, query):
        """
        Searches sample sequence for query subsequence
        :param query: Dict containing 'sequence' DNA string entry
        :return: list of mappy alignment objects
        """
        sequence = query['sequence']
        alignments = list(self.aligner.map(sequence))
        chrom_alignments = [a for a in alignments if a.ctg.startswith(query['chrom'])]

        return chrom_alignments

    def score_all(self, location_tolerance=float('inf'), phred_threshold=50):
        """
        For each SV, checks whether a subsequence matching its result exists in the sample sequence
        Breaks down accuracy by SV type and total
        :param location_tolerance: maximum base pair location distance to count as a correct match (default infinite)
        :param phred_threshold: maximum phred score to count as a correct match
        :return: rate of correct matches (number matches / total SVs), i.e., precision
        """

        total_calls = Counter()
        correct_calls = Counter()
        precision = {}

        for svid in self.variants.keys():
            sv_type = self.variants[svid][0].info['SVTYPE']
            total_calls[sv_type] += 1

            sequences = self.simulate_subsequences(svid)
            matched = []

            # check that all subsequences match
            for sequence in sequences:
                alignments = self.match_subsequence(sequence)

                close_alignments = [a for a in alignments if abs(a.r_st - sequence['location']) <= location_tolerance]
                matched.append(any([a.mapq > phred_threshold for a in close_alignments]))

            if all(matched) and len(matched) == len(sequences):
                correct_calls[sv_type] += 1

            precision[sv_type] = correct_calls[sv_type] / total_calls[sv_type]


        precision['ALL'] = sum(correct_calls.values()) / sum(total_calls.values())

        return precision

def main():
    parser = argparse.ArgumentParser(description='Score VCF SV calls against a reference and sample genome')
    parser.add_argument('--reference', help='Reference genome .fa file', dest='reference')
    parser.add_argument('--sample', help='Sample genome .fa file', dest='sample')
    parser.add_argument('--calls', help='VCF containing called SVs', dest='calls')
    parser.add_argument('--output_dir', help='Output directory', dest='output_dir')
    parser.add_argument('--buffer', help='Subsequence context buffer', type=int, dest='buffer', default=500)
    args = parser.parse_args()

    scorer = AlignScorer(args.reference, args.sample, args.calls, args.output_dir, args.buffer)
    sequences = scorer.simulate_subsequences('sv0')
    scorer.match_subsequence(sequences[0])

    precision = scorer.score_all(phred_threshold=50)

    print(f"Precision: {precision}")


if __name__ == '__main__':
    main()

'''
--reference /data/refs/refdata-GRCh38-2.1.0/fasta/genome.fa --sample /data/refs/HG002/hg002v1.1.fasta --calls ./test_calls.vcf
--reference /data/refs/refdata-GRCh38-2.1.0/fasta/genome.chr21.fa --sample /data/bert/groovi-infra/sim_data/sim.hapA.fa --calls /data/bert/groovi-infra/sim_data/sim.vcf --output_dir output --buffer 50
'''