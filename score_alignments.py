import argparse
import mappy
import pysam
import os
import tempfile
from intervaltree import IntervalTree
from itertools import groupby
from tqdm import tqdm

from collections import defaultdict, Counter
from typing import Dict, List, Tuple

from pysam.libcvcf import VCFRecord

import logging

logger = logging.getLogger(__name__)

def get_start_stop(rec: VCFRecord) -> Tuple[int, int]:
    """
    Converts VCF start/stop coordinates to Python style
    VCF range coordinates are 1-indexed and have inclusive ends
    :param rec: VCF record with start and end
    :return: Tuple containing start and stop coordinates
    """
    start = rec.start - 1
    stop = rec.stop

    return start, stop

class AlignScorer(object):
    def __init__(self, ref_sequence, sample_sequence, callset_vcf, output_dir, buffer, gap_file):
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

        self.blacklist = defaultdict(IntervalTree)

        if gap_file:
            self.load_blacklist(gap_file)

        print(f"Loading reference into aligner")
        self.aligner = mappy.Aligner(fn_idx_in=sample_sequence, preset='map-pb')
        print(f"Loaded aligner with sequences {self.aligner.seq_names}")

        print("Loading sample")
        self.ref_fasta = pysam.FastaFile(ref_sequence)
        print("Sample loaded")

        self.sample_sequence = sample_sequence
        self.working_dir = output_dir
        self.buffer = buffer

    def load_blacklist(self, gap_file):
        with open(gap_file, 'r') as f:
            for line in f:
                row = line.strip().split()

                chrom = row[1]
                start, stop = int(row[2]), int(row[3])
                region_type = row[7]

                print(f"{chrom}\t{start}\t{stop}\t{region_type}")

                if region_type in ['telomere', 'centromere']:
                    self.blacklist[chrom][start:stop] = region_type

        print("Loaded blacklist")

        # TODO: for now, leaving this unused. I will reimplement this in the stitching process. Eventually should include here independently as well.


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


        # Reorder records to place any insertions at the end, sorted by target location
        target_records = []
        in_place_records = []

        for rec in records:
            if 'TARGET' in rec.info:
                target_records.append(rec)
            else:
                in_place_records.append(rec)

        if len(target_records) > 1:
            target_records = sorted(target_records, key=lambda rec: (rec.info['TARGET'], rec.info.get('INSORD', 0)), reverse=True)
            # print(f"WARNING: The list for {svid}-{sv_type} contains more than one record with target. Sorting in reverse order by target and INSORD if available, but conflicts are possible")

        records = in_place_records + target_records

        offset = max(0, min([rec.start for rec in records] + [rec.info['TARGET'] for rec in target_records]) - self.buffer)
        sequence_end = max([rec.stop for rec in records] + [rec.info['TARGET'] for rec in target_records]) + self.buffer

        orig_sequence = list(self.ref_fasta.fetch(reference=chrom)[offset:sequence_end])
        new_sequence = orig_sequence.copy()
        changed_mask = [False] * len(new_sequence)

        delete_placeholder = ''
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', '': '',
                      'a': 't', 'c': 'g', 'g': 'c', 't': 'a'}  # preserve soft masking

        def reverse_complement(seq: List):
            return [complement[x] for x in reversed(seq)]


        queries = []

        for rec in records:
            # Supported operations:
            # CUT (aka DEL)
            # COPY-PASTE (aka dDUP, DUP)
            # CUT-PASTE (aka nrTRA)
            # COPYinv_PASTE (aka INV_dDUP)
            # CUTinv-PASTE (aka INV-nrTRA)

            start, stop = get_start_stop(rec)
            target = rec.info.get('TARGET', stop + 1) - offset  # Convert target to 0-index
            start -= offset
            stop -= offset

            if rec.info.get('TARGET_CHROM', chrom) != chrom:
                print(f"WARNING: Skipping record: {rec.id}-{sv_type} because interchromosome target. Interchromosome checks not implemented yet.")
                return []


            if rec.info['OP_TYPE'] == 'CUT' or rec.info['SVTYPE'] == 'DEL':
                new_sequence[start:stop] = [delete_placeholder] * (stop - start)

                changed_mask[start:stop] = [True] * (stop - start)
            elif rec.info['OP_TYPE'] == 'INV' or rec.info['SVTYPE'] == 'INV':
                new_sequence[start:stop] = reverse_complement(orig_sequence[start:stop])

                changed_mask[start:stop] = [True] * (stop - start)
            elif rec.info['OP_TYPE'] == 'COPY-PASTE' or rec.info['SVTYPE'] in ['DUP', 'dDUP']:
                clip = orig_sequence[start:stop]
                new_sequence = new_sequence[:target] + clip + new_sequence[target:]

                changed_mask = changed_mask[:target] + [True] * len(clip) + changed_mask[target:]
            elif rec.info['OP_TYPE'] in ['CUT-PASTE'] or rec.info['SVTYPE'] == 'nrTRA':
                clip = orig_sequence[start:stop]
                new_sequence[start:stop] = [delete_placeholder] * (stop - start)
                changed_mask[start:stop] = [True] * len(clip)

                new_sequence = new_sequence[:target] + clip + new_sequence[target:]
                changed_mask = changed_mask[:target] + [True] * len(clip) + changed_mask[target:]
            elif rec.info['OP_TYPE'] == 'COPYinv-PASTE' or rec.info['SVTYPE'] in ['INV_dDUP']:
                clip = reverse_complement(orig_sequence[start:stop])
                new_sequence = new_sequence[:target] + clip + new_sequence[target:]

                changed_mask = changed_mask[:target] + [True] * len(clip) + changed_mask[target:]
            elif rec.info['OP_TYPE'] == 'CUTinv-PASTE' or rec.info['SVTYPE'] in ['INV_nrTRA']:
                clip = reverse_complement(orig_sequence[start:stop])

                new_sequence[start:stop] = [delete_placeholder] * (stop - start)
                new_sequence = new_sequence[:target] + clip + new_sequence[target:]

                changed_mask[start:stop] = [True] * len(clip)
                changed_mask = changed_mask[:target] + [True] * len(clip) + changed_mask[target:]
            else:
                print(f"WARNING: Unknown OP_TYPE: {rec.info['OP_TYPE']}")

        if len(changed_mask) == 0:
            print(f"WARNING: No changed intervals found for {svid}-{sv_type}")
            return []

        intervals = []

        # Get contiguous blocks of changed indices in the subsequence
        i = 0
        for value, group in groupby(changed_mask):
            group_len = len(list(group))
            if value:
                intervals.append((i, i + group_len))
            i += group_len

        for start, stop in intervals:
            adjusted_start = max(0, start - self.buffer)
            adjusted_stop = min(stop + self.buffer + 1, len(new_sequence))
            sequence = ''.join(new_sequence[adjusted_start:adjusted_stop])
            query = {
                'chrom': rec.chrom,
                'svtype': sv_type,
                'svid': svid,
                'sequence': sequence,
                'location': adjusted_start + offset,
                'length': len(sequence),
            }

            queries.append(query)

        # print(f'Checking {svid}-{sv_type} with {len(queries)} queries')

        return queries

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


    def score_all(self, location_tolerance=float('inf'), error_threshhold=0.1):
        """
        For each SV, checks whether a subsequence matching its result exists in the sample sequence
        Breaks down accuracy by SV type and total
        :param location_tolerance: maximum base pair location distance to count as a correct match (default infinite)
        :param error_threshold: maximum fraction of string mismatch
        :return: rate of correct matches (number matches / total SVs), i.e., precision
        """

        total_calls = Counter()
        correct_calls = Counter()
        precision = {}

        def check_match(alignment, query):
            """
            Calculate the independent error rate of a mappy alignment with the query sequence (ignoring redundancy penalties)
            :param alignment: mappy Alignment object
            :param query: string containing query sequence
            :return: Boolean indicating whether error rate is less than error_threshold
            """
            aligned_query_segment_length = alignment.q_en - alignment.q_st
            len_unaligned = len(query) - aligned_query_segment_length

            # Numerator: Internal errors (NM) + Penalty for unaligned ends
            nm_total = alignment.NM + len_unaligned

            # Denominator: Aligned block length (blen) + Length of unaligned ends
            len_norm = alignment.blen + len_unaligned

            error_rate = nm_total / len_norm if len_norm > 0 else 1.0

            return error_rate <= error_threshhold

        overall_count = 0
        overall_correct = 0

        # # Filter for debugging
        # debug_examples = [
        #     'sv8',
        # ]
        # self.variants = {key:self.variants[key] for key in debug_examples}

        pbar = tqdm(self.variants.keys(), total=len(self.variants), desc=f'Scoring SVs')

        for svid in pbar:
            sv_type = self.variants[svid][0].info['SVTYPE']
            total_calls[sv_type] += 1

            sequences = self.simulate_subsequences(svid)
            matched = []

            # check that all subsequences match
            for sequence in sequences:
                alignments = self.match_subsequence(sequence)

                close_alignments = [a for a in alignments if abs(a.r_st - sequence['location']) <= location_tolerance]
                matched.append(any([check_match(a, sequence['sequence']) for a in close_alignments]))
                # print(f"Query match for {sequence['svid']}-{sequence['svtype']}: {matched[-1]}")

            coords = sequences[0]['location'], sequences[0]['location'] + sequences[0]['length']

            if len(sequences) > 0 and all(matched) and len(matched) == len(sequences):
                # print("Match successful")
                correct_calls[sv_type] += 1
                overall_correct += 1
                hit_miss = 'hit'
            else:
                hit_miss = 'miss'
                # print("Match not found")
            logger.info(f'{svid}\t{sv_type}\t{sequences[0]["chrom"]}:{coords[0]}-{coords[1]}\t{hit_miss}\n')

            overall_count += 1

            pbar.set_description(f'Scoring SVs. Current precision {overall_correct / overall_count:.2f} ({overall_correct} / {overall_count})')

            precision[sv_type] = correct_calls[sv_type] / total_calls[sv_type]


        precision['ALL'] = sum(correct_calls.values()) / sum(total_calls.values())
        correct_calls['ALL'] = sum(correct_calls.values())
        total_calls['ALL'] = sum(total_calls.values())

        return precision, correct_calls, total_calls

def main():
    logging.basicConfig(filename='log.txt', level=logging.INFO, filemode='w')

    parser = argparse.ArgumentParser(description='Score VCF SV calls against a reference and sample genome')
    parser.add_argument('--reference', help='Reference genome .fa file', dest='reference')
    parser.add_argument('--sample', help='Sample genome .fa file', dest='sample')
    parser.add_argument('--calls', help='VCF containing called SVs', dest='calls')
    parser.add_argument('--output_dir', help='Output directory', dest='output_dir')
    parser.add_argument('--buffer', help='Subsequence context buffer', type=int, dest='buffer', default=500)
    parser.add_argument('--gap_file', help='Tab-delimited file containing centromere and telomere regions', default=None)
    args = parser.parse_args()

    scorer = AlignScorer(args.reference, args.sample, args.calls, args.output_dir, args.buffer, args.gap_file)
    # sequences = scorer.simulate_subsequences('sv0')
    # scorer.match_subsequence(sequences[0])

    precision, correct_calls, total_calls = scorer.score_all()

    print(f"Precision: {precision}")

    import pandas as pd

    df = pd.DataFrame({
        'correct_calls': correct_calls,
        'total_calls': total_calls,
        'precision': precision,
    })

    print(df)




if __name__ == '__main__':
    main()

'''
cat ./sim_data/sim.hapA.fa ./sim_data/sim.hapB.fa > ./sim_data/sim.combined.fa

--reference ./data/genome.chr21.fa --sample ./sim_data/sim.combined.fa --calls ./sim_data/sim.vcf --output_dir output --buffer 500

real data settings
--reference /Users/huangber/remote/data/refs/refdata-GRCh38-2.1.0/fasta/genome.fa --sample /Users/huangber/remote/data/refs/HG002/hg002v1.1.fasta --calls ../groovi/output.vcf --output_dir output --buffer 500
--reference ./data/hg19.genome.fa --sample ./data/hg002v1.1.fasta --calls ../groovi/data/output.vcf --output_dir output --buffer 500
--reference /data/refs/refdata-GRCh38-2.1.0/fasta/genome.fa --sample /data/refs/HG002/hg002v1.1.fasta --calls ../groovi/output.vcf --output_dir output --buffer 500
--reference ./data/hg19.genome.fa --sample ./data/hg002.mmi --calls ../groovi/data/output.vcf --output_dir output --buffer 1000
'''

# TODO:
#  - turn on filtering,
#  - set up whitelist/blacklist for centromeres/telomeres,
#  - check reference for pseudo-precision
#  - make sure merged classes are properly handled
#  -