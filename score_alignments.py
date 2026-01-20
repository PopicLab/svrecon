import argparse
import mappy
import pysam
import os
import tempfile
from intervaltree import IntervalTree
from itertools import groupby

import traceback
from pysam.libcfaidx import FastaFile
from tqdm import tqdm
import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading
import yaml
from pathlib import Path

from collections import defaultdict, Counter
from typing import Dict, List, Tuple

from pysam import VariantRecord

import logging

logger = logging.getLogger(__name__)

def get_start_stop(rec: VariantRecord) -> Tuple[int, int]:
    """
    Converts VCF start/stop coordinates to Python style
    VCF range coordinates are 1-indexed and have inclusive ends
    :param rec: VCF record with start and end
    :return: Tuple containing start and stop coordinates
    """
    start = rec.start - 1
    stop = rec.stop

    return start, stop

# Global variable to hold the sample Aligner and reference FastaFile instance in each worker process
# (Must be accessible by the worker function)
global_aligner = None
global_ref = None


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

    return error_rate

def simulate_subsequences(records: List[VariantRecord], buffer: int) -> List[Dict]:
    """
    Execute operations defined by called SV on the reference genome
    :param: records: list of VCF records describing operations
    :param: buffer: buffer size for surrounding context
    :param: ref_fasta: Fasta file containing reference genome
    :return: List of sequence dictionaries, each element containing "chrom", "sequence", and "location"

    For SVs without dispersions, the list should only contain a single element. For dispersions,
    the list should include the evidence of the SV at the source and the target.
    """
    global global_ref

    sv_type = records[0].info['SVTYPE']
    svid = records[0].info['SVID']
    chrom = records[0].chrom

    # Reorder records to place any insertions at the end, sorted by target location
    target_records = []
    in_place_records = []

    for rec in records:
        # Patch to fill in DUP targets
        if rec.info['SVTYPE'] == 'DUP':
            rec.info['TARGET'] = rec.stop

        if 'TARGET' in rec.info:
            target_records.append(rec)
        else:
            in_place_records.append(rec)

    if len(target_records) > 1:
        target_records = sorted(target_records, key=lambda rec: (rec.info['TARGET'], rec.info.get('INSORD', 0)), reverse=True)
        # print(f"WARNING: The list for {svid}-{sv_type} contains more than one record with target. Sorting in reverse order by target and INSORD if available, but conflicts are possible")

    records = in_place_records + target_records

    offset = max(0, min([rec.start for rec in records]) - buffer)
    sequence_end = max([rec.stop for rec in records]) + buffer

    target_offset = max(0, min([rec.info['TARGET'] for rec in target_records]) - buffer) if target_records else 0
    target_sequence_end = max([rec.info['TARGET'] for rec in target_records]) + buffer if target_records else 0

    merged_target_sequence = offset <= target_offset + buffer and target_sequence_end - buffer <= sequence_end
    if merged_target_sequence:
        target_offset = offset

    orig_sequence = list(global_ref[chrom][offset:sequence_end].decode('ascii'))
    new_sequence = orig_sequence.copy()
    changed_mask = [False] * len(new_sequence)

    if not merged_target_sequence:
        new_target_sequence = list(global_ref[chrom][target_offset:target_sequence_end].decode('ascii'))
        target_changed_mask = [False] * len(new_target_sequence)

    delete_placeholder = ''
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', '': '',
                  'a': 't', 'c': 'g', 'g': 'c', 't': 'a',   # preserve soft masking
                  'N': 'N', 'n': 'n'}

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
        target = rec.info.get('TARGET', stop + 1) - target_offset  # Convert target to 0-index and offset
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
            if merged_target_sequence:
                new_sequence = new_sequence[:target] + clip + new_sequence[target:]
                changed_mask = changed_mask[:target] + [True] * len(clip) + changed_mask[target:]
            else:
                new_target_sequence = new_target_sequence[:target] + clip + new_target_sequence[target:]
                target_changed_mask = target_changed_mask[:target] + [True] * len(clip) + target_changed_mask[target:]

        elif rec.info['OP_TYPE'] in ['CUT-PASTE'] or rec.info['SVTYPE'] == 'nrTRA':
            clip = orig_sequence[start:stop]
            new_sequence[start:stop] = [delete_placeholder] * (stop - start)
            changed_mask[start:stop] = [True] * len(clip)

            if merged_target_sequence:
                new_sequence = new_sequence[:target] + clip + new_sequence[target:]
                changed_mask = changed_mask[:target] + [True] * len(clip) + changed_mask[target:]
            else:
                new_target_sequence = new_target_sequence[:target] + clip + new_target_sequence[target:]
                target_changed_mask = target_changed_mask[:target] + [True] * len(clip) + target_changed_mask[target:]
        elif rec.info['OP_TYPE'] == 'COPYinv-PASTE' or rec.info['SVTYPE'] in ['INV_dDUP']:
            clip = reverse_complement(orig_sequence[start:stop])

            if merged_target_sequence:
                new_sequence = new_sequence[:target] + clip + new_sequence[target:]
                changed_mask = changed_mask[:target] + [True] * len(clip) + changed_mask[target:]
            else:
                new_target_sequence = new_target_sequence[:target] + clip + new_target_sequence[target:]
                target_changed_mask = target_changed_mask[:target] + [True] * len(clip) + target_changed_mask[target:]
        elif rec.info['OP_TYPE'] == 'CUTinv-PASTE' or rec.info['SVTYPE'] in ['INV_nrTRA']:
            clip = reverse_complement(orig_sequence[start:stop])

            new_sequence[start:stop] = [delete_placeholder] * (stop - start)
            changed_mask[start:stop] = [True] * len(clip)

            if merged_target_sequence:
                new_sequence = new_sequence[:target] + clip + new_sequence[target:]
                changed_mask = changed_mask[:target] + [True] * len(clip) + changed_mask[target:]
            else:
                new_target_sequence = new_target_sequence[:target] + clip + new_target_sequence[target:]
                target_changed_mask = target_changed_mask[:target] + [True] * len(clip) + target_changed_mask[target:]
        else:
            print(f"WARNING: Unknown OP_TYPE: {rec.info['OP_TYPE']}")

    if len(changed_mask) == 0:
        print(f"WARNING: No changed intervals found for {svid}-{sv_type}")
        return []

    MERGE_TOLERANCE = 5  # merge intervals that are within this many base pairs (usually alignment or boundary case error)

    queries.extend(get_changed_subsequences(new_sequence, changed_mask, MERGE_TOLERANCE, offset, buffer, svid, chrom, sv_type))
    if not merged_target_sequence:
        queries.extend(get_changed_subsequences(new_target_sequence, target_changed_mask, MERGE_TOLERANCE, target_offset, buffer, svid, chrom, sv_type))

    # print(f'Checking {svid}-{sv_type} with {len(queries)} queries')

    return queries


def get_changed_subsequences(new_sequence, changed_mask, tolerance, offset, buffer, svid, chrom, sv_type):
    changed_intervals = []

    queries = []

    # Get contiguous blocks of changed indices in the subsequence
    current_index = 0
    for value, group in groupby(changed_mask):
        group_len = len(list(group))
        if value:
            interval_start = current_index
            prev_interval = changed_intervals[-1] if changed_intervals else None

            if prev_interval and current_index - prev_interval[1] <= tolerance:
                prev_interval = changed_intervals.pop()
                interval_start = prev_interval[0]

            changed_intervals.append((interval_start, current_index + group_len))
        current_index += group_len


    for start, stop in changed_intervals:
        adjusted_start = max(0, start - buffer)
        adjusted_stop = min(stop + buffer + 1, len(new_sequence))
        sequence = ''.join(new_sequence[adjusted_start:adjusted_stop])
        query = {
            'chrom': chrom,
            'svtype': sv_type,
            'svid': svid,
            'sequence': sequence,
            'location': adjusted_start + offset,
            'length': len(sequence),
        }
        queries.append(query)

    return queries




def score_sv(records: List[VariantRecord], buffer: int, location_tolerance: int, error_threshold=0.1):
    """

    :param records:
    :param buffer:
    :param location_tolerance:
    :param error_threshold:
    :return:
    """
    try:
        svid = records[0].info['SVID']
        sv_type = records[0].info['SVTYPE']

        sequences = simulate_subsequences(records, buffer)
        match_scores = []

        # check that all subsequences match
        for query in sequences:
            sequence = query['sequence']
            all_alignments = list(global_aligner.map(sequence))

            alignments = [a for a in all_alignments if a.ctg.startswith(query['chrom'])]

            alignments = [a for a in alignments if abs(a.r_st - query['location']) <= location_tolerance]

            if len(alignments) > 0:
                match_scores.append(min([check_match(a, sequence) for a in alignments]))
            else:
                match_scores.append(1)  # no candidate alignments found, marking this as a miss by adding a match error of 1.0 (maximum error)

        coords = sequences[0]['location'], sequences[0]['location'] + sequences[0]['length']

        result = {}

        if len(sequences) > 0 and all([score <= error_threshold for score in match_scores]):
            # print("Match successful")
            is_correct = True
            hit_miss = 'hit'
        else:
            hit_miss = 'miss'
            is_correct = False
            # print("Match not found")

        return {
            'svid': svid,
            'sv_type': sv_type,
            'is_correct': is_correct,
            'coords': coords,
            'match_scores': match_scores,
            'sequences': sequences,
        }
    except Exception as e:
        # 1. Log the full traceback *immediately*
        error_msg = f"Worker failed for SVID: {records[0].info.get('SVID', 'Unknown')}. Error: {e}\n{traceback.format_exc()}"
        print(error_msg)

        # 2. Re-raise the exception so the main process still knows it failed
        raise e


class AlignScorer(object):
    def __init__(self, callset_vcf, output_dir, buffer, gap_file):
        """
        :param callset_vcf: File containing SV callset
        :param output_dir: Output directory for auxiliary files
        :param buffer: Context buffer length
        """
        self.variants = None
        self.vcf_header = None
        self.read_vcf(callset_vcf)

        if gap_file:
            self.exclude_list = self.load_exclude_list(gap_file)

            filtered_variants = []
            for svid, records in self.variants.items():
                allowed = True
                for rec in records:
                    target_chrom = rec.info['TARGET_CHROM'] if 'TARGET_CHROM' in rec.info else None
                    if self.exclude_list[rec.chrom].overlap(rec.start, rec.stop) or \
                        target_chrom and self.exclude_list[target_chrom].overlap(rec.info['TARGET'], rec.info['TARGET'] + 1):
                        # record is in the excluded regions
                        allowed = False
                        break
                if allowed:
                    filtered_variants.append((svid, records))

            self.variants = {key: value for (key, value) in filtered_variants}

        self.working_dir = output_dir
        self.buffer = buffer

    def load_exclude_list(self, gap_file):
        exclude_list = defaultdict(IntervalTree)
        with open(gap_file, 'r') as f:
            for line in f:
                row = line.strip().split()

                chrom = row[1]
                start, stop = int(row[2]), int(row[3])
                region_type = row[7]

                print(f"{chrom}\t{start}\t{stop}\t{region_type}")

                # if region_type in ['telomere', 'centromere']:
                # Include all regions in this gap file
                exclude_list[chrom][start:stop] = region_type

        return exclude_list


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


    def score_all(self, location_tolerance=float('inf'), error_threshold=0.1, n_threads=16):
        """
        For each SV, checks whether a subsequence matching its result exists in the sample sequence
        Breaks down accuracy by SV type and total
        :param location_tolerance: maximum base pair location distance to count as a correct match (default infinite)
        :param error_threshold: maximum fraction of string mismatch
        :return: rate of correct matches (number matches / total SVs), i.e., precision
        """
        total_calls = Counter()
        correct_calls = Counter()
        overall_count = 0
        overall_correct = 0
        precision = {}
        results_list = []

        with ThreadPoolExecutor(max_workers=n_threads) as executor:
            # Submit all SV tasks to the pool
            futures = {
                executor.submit(
                    score_sv,
                    self.variants[svid],
                    self.buffer,
                    location_tolerance,
                    error_threshold
                )
                for svid in self.variants.keys()
            }

            # Aggregate results and manage the progress bar
            pbar = tqdm(as_completed(futures), total=len(self.variants), desc=f'Scoring SVs')

            for future in pbar:
                try:
                    result = future.result()

                    sv_type = result['sv_type']
                    svid = result['svid']
                    sequences = result['sequences']
                    coords = result['coords']
                    hit_miss = 'hit' if result['is_correct'] else 'miss'
                    match_scores = result['match_scores']
                    logger.info(f'{svid}\t{sv_type}\t{sequences[0]["chrom"]}:{coords[0]}-{coords[1]}\t{hit_miss}\t{match_scores}')

                    total_calls[sv_type] += 1
                    overall_count += 1

                    if result['is_correct']:
                        correct_calls[sv_type] += 1
                        overall_correct += 1
                    else:
                        correct_calls[sv_type] += 0

                    # Update display
                    pbar.set_description(
                        f'Scoring SVs. Current precision {overall_correct / overall_count:.2f} ({overall_correct} / {overall_count})')
                    precision[sv_type] = correct_calls[sv_type] / total_calls[sv_type]

                except Exception as exc:
                    logger.error(f'SV processing generated an exception: {exc}')

        # 4. Final aggregation
        precision['ALL'] = sum(correct_calls.values()) / sum(total_calls.values())
        correct_calls['ALL'] = sum(correct_calls.values())
        total_calls['ALL'] = sum(total_calls.values())

        return precision, correct_calls, total_calls

def export_igv_session(calls, bam, classified, timestamp):
    # For now, hard coded to hg19 as reference
    import xml.etree.ElementTree as ET

    output_filename = f"./igv_sessions/session_{timestamp}.xml"

    # 2. Build XML
    session = ET.Element("Session", genome="hg19", version="8")
    resources = ET.SubElement(session, "Resources")

    ET.SubElement(resources, "Resource", path=calls, type="vcf")
    ET.SubElement(resources, "Resource", path=classified, type="vcf")
    ET.SubElement(resources, "Resource", path=bam, type="bam")

    # 3. Write to file with standard header
    tree = ET.ElementTree(session)
    tree.write(output_filename, encoding="utf-8", xml_declaration=True)

    print(f"File saved to {output_filename}")

def main():
    timestamp = datetime.datetime.now().strftime('%Y-%m-%d_%H%M%S')
    log_filename = os.path.join("./logs", f'reconstruction_{timestamp}.log')
    logging.basicConfig(filename=log_filename, level=logging.DEBUG, filemode='w')

    parser = argparse.ArgumentParser(description='Score VCF SV calls against a reference and sample genome')
    parser.add_argument('--reference', help='Reference genome .fa file', dest='reference')
    parser.add_argument('--sample', help='Sample genome .fa file', dest='sample')
    parser.add_argument('--calls', help='VCF containing called SVs', dest='calls')
    parser.add_argument('--output_dir', help='Output directory', dest='output_dir')
    parser.add_argument('--buffer', help='Subsequence context buffer', type=int, dest='buffer', default=500)
    parser.add_argument('--gap_file', help='Tab-delimited file containing centromere and telomere regions', default=None)
    parser.add_argument('--config', help='Groovi call config used to infer other params', dest='config')
    args = parser.parse_args()

    logger.info(f'Config: {vars(args)}')

    if args.config:
        logger.info(f'Inferring params from {args.config}')
        with open(args.config, 'r') as file:
            config_data = yaml.safe_load(file)

            # Find groovi output
            experiment_dir = str(Path(args.config).parent.resolve())
            results_dir = os.path.join(experiment_dir, "results")
            if args.calls is None:
                args.calls = os.path.join(results_dir, 'groovi.vcf')

            # Attempt to find InsilicoSV sample genome
            # Assumes a common path organization of mounted drives.
            exp_path = Path(experiment_dir)
            prefix = data_dir = next(p for p in exp_path.parents if p.name == 'data').parent
            bam_path = Path(config_data['bam'])
            data_index = bam_path.parts.index('data')
            bam_path = Path(*bam_path.parts[data_index:])

            bam = prefix / bam_path

            args.bam = str(bam)


            sample = bam.parent.parent / 'VCF/sim.fa'

            if sample.exists() and args.sample is None:
                args.sample = str(sample)

            # Grab reference from config
            fa_path = Path(config_data['fa'])
            data_index = fa_path.parts.index('data')
            fa_path = Path(*fa_path.parts[data_index:])

            args.classified = os.path.join(experiment_dir, "results/groovi_bkps_classified.vcf")

            if args.reference is None:
                args.reference = str(prefix / fa_path)

            logger.info(f'Config: {vars(args)}')


    global global_ref
    global global_aligner

    print('Initializing scorer and loading callset')
    scorer = AlignScorer(args.calls, args.output_dir, args.buffer, args.gap_file)

    print('Finding relevant chromosomes')
    chroms = set()
    for records in scorer.variants.values():
        for record in records:
            chroms.add(record.chrom)
            if 'TARGET_CHROM' in record.info:
                chroms.add(record.info['TARGET_CHROM'])
    print(f'Found {len(chroms)} referenced chromosomes in callset')

    print("Loading reference")
    ref_data = {}
    with pysam.FastaFile(args.reference) as f:
        for chrom in chroms:
            sequence_string = f.fetch(reference=chrom)
            # ref_data[chrom] = list(sequence_string)
            ref_data[chrom] = bytearray(sequence_string, 'ascii')
    global_ref = ref_data
    print(f"Reference pre-loaded with {len(global_ref)} chromosomes.")

    print("Loading index into aligner")
    global_aligner = mappy.Aligner(args.sample, preset='map-pb', n_threads=16)
    print("Aligner ready")

    # # Filter for debugging
    # debug_examples = [
    #     'sv8',
    #     'sv384',
    # ]
    # scorer.variants = {key:scorer.variants[key] for key in debug_examples}
    # scorer.variants = {key:scorer.variants[key] for key in scorer.variants if scorer.variants[key][0].info['SVTYPE'] in ['dupINVdup']}

    precision, correct_calls, total_calls = scorer.score_all()

    import pandas as pd

    df = pd.DataFrame({
        'correct_calls': correct_calls,
        'total_calls': total_calls,
        'precision': precision,
    })

    print(df)
    logger.info(f'Score table: {df}')

    export_igv_session(args.calls, args.bam, args.classified, timestamp)


if __name__ == '__main__':
    main()

'''
cat ./sim_data/sim.hapA.fa ./sim_data/sim.hapB.fa > ./sim_data/sim.combined.fa

--reference ./data/genome.chr21.fa --sample ./sim_data/sim.combined.fa --calls ./sim_data/sim.vcf --output_dir output --buffer 500

# sim benchmark tests
--reference ./data/genome.fa.chr1.fa --sample ./sim_data/long_sim.fa --calls ./sim_data/long_sim.vcf --output_dir output --buffer 500


real data settings
--reference /Users/huangber/remote/data/refs/refdata-GRCh38-2.1.0/fasta/genome.fa --sample /Users/huangber/remote/data/refs/HG002/hg002v1.1.fasta --calls ../groovi/output.vcf --output_dir output --buffer 500
--reference ./data/hg19.genome.fa --sample ./data/hg002v1.1.fasta --calls ../groovi/data/output.vcf --output_dir output --buffer 500
--reference /data/refs/refdata-GRCh38-2.1.0/fasta/genome.fa --sample /data/refs/HG002/hg002v1.1.fasta --calls ../groovi/output.vcf --output_dir output --buffer 500
--reference ./data/hg19.genome.fa --sample ./data/hg002.mmi --calls ../groovi/data/output.vcf --output_dir output --buffer 1000
'''