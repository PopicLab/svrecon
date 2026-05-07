import argparse
import datetime
import hashlib
import logging
import os
import sys
import tempfile
import traceback
from collections import defaultdict, Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from itertools import groupby
from pathlib import Path
from typing import Dict, List, Tuple, Union

import edlib
import mappy
import pysam
import yaml
from intervaltree import IntervalTree
from pysam import VariantRecord
from tqdm import tqdm

logger = logging.getLogger(__name__)

MIN_EDLIB_QUERY = 10000

def get_chrom_aligner(sample_fasta: str, chrom: str, cache_dir: str, align_params: dict,
                      threads: int = 4) -> mappy.Aligner:
    """Loads a per-chromosome MMI index, building it first if it doesn't exist."""
    mmi_path = os.path.join(cache_dir, f"{chrom}.mmi")

    # If we already have the index in the cache, load it directly
    if os.path.exists(mmi_path):
        logger.info(f"Loading existing index for {chrom} from {mmi_path}")
        return mappy.Aligner(mmi_path, **align_params)

    # Build the index by extracting the relevant contigs
    logger.info(f"Index not found for {chrom}. Extracting sequences...")
    chrom_fa_path = os.path.join(cache_dir, f"{chrom}.fa")
    written_contigs = 0

    with pysam.FastaFile(sample_fasta) as fasta:
        with open(chrom_fa_path, 'w') as out_f:
            for ref in fasta.references:
                # Catch exact chrom or haplotype extensions (e.g., chr1, chr1_paternal)
                if ref == chrom or ref.startswith(f"{chrom}_"):
                    seq = fasta.fetch(ref)
                    out_f.write(f">{ref}\n{seq}\n")
                    written_contigs += 1

    if written_contigs == 0:
        if os.path.exists(chrom_fa_path):
            os.remove(chrom_fa_path)
        return None

    # Build the mappy index, saving it to mmi_path
    logger.info(f"Building MMI index for {chrom} at {mmi_path}...")
    aligner = mappy.Aligner(chrom_fa_path, n_threads=threads, fn_idx_out=mmi_path, **align_params)

    # Clean up the intermediate FASTA to save disk space
    if os.path.exists(chrom_fa_path):
        os.remove(chrom_fa_path)

    return aligner


def update_args_from_config(args):
    """
    Infer parameters from groovi config files
    :param args: args passed into this tool including args.config pointing to a groovi config file
    :return: updated args object
    """
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
        if not bam_path:
            data_index = bam_path.parts.index('data')
            bam_path = Path(*bam_path.parts[data_index:])

        bam = prefix / bam_path

        if not args.bam:
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

        logger.info(f'Config updated: {vars(args)}')
    return args


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


def reverse_complement(seq: Union[str, List[str]]) -> List[str]:
    """Returns the reverse complement of a given DNA sequence string."""
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
                  'a': 't', 'c': 'g', 'g': 'c', 't': 'a',
                  'N': 'N', 'n': 'n'}
    return [complement[x] for x in reversed(seq)]


# Global variable to share the sample and reference for each worker process
# (Must be accessible by the worker function)
global_ref = None
global_sample = None
global_aligners = None


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
        target_records = sorted(target_records, key=lambda rec: (rec.info['TARGET'], rec.info.get('INSORD', 0)),
                                reverse=True)

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

    queries = []

    for rec in records:
        start, stop = get_start_stop(rec)
        target = rec.info.get('TARGET', stop + 1) - target_offset  # Convert target to 0-index and offset
        start -= offset
        stop -= offset

        if rec.info.get('TARGET_CHROM', chrom) != chrom:
            logger.warning(
                f"Skipping record: {rec.id}-{sv_type} because interchromosome target. Interchromosome checks not implemented yet.")
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
            logger.warning(f"Unknown OP_TYPE: {rec.info['OP_TYPE']}")

    if len(changed_mask) == 0:
        logger.warning(f"No changed intervals found for {svid}-{sv_type}")
        return []

    MERGE_TOLERANCE = 5  # merge intervals that are within this many base pairs (usually alignment or boundary case error)

    queries.extend(
        get_changed_subsequences(new_sequence, changed_mask, MERGE_TOLERANCE, offset, buffer, svid, chrom, sv_type))
    if not merged_target_sequence:
        queries.extend(
            get_changed_subsequences(new_target_sequence, target_changed_mask, MERGE_TOLERANCE, target_offset, buffer,
                                     svid, chrom, sv_type))

    return queries


def get_changed_subsequences(new_sequence, changed_mask, tolerance, offset, buffer, svid, chrom, sv_type):
    """
    :param new_sequence:
    :param changed_mask:
    :param tolerance:
    :param offset:
    :param buffer:
    :param svid:
    :param chrom:
    :param sv_type:
    :return:
    """
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


def load_fasta_to_bytes(filename: str, chroms) -> Dict[str, bytearray]:
    data = {}
    with pysam.FastaFile(filename) as f:
        fasta_refs = f.references
        for chrom in chroms:
            chrom_lower = chrom.lower()
            for ref in fasta_refs:
                ref_lower = ref.lower()
                if ref_lower == chrom_lower or ref_lower.startswith(f"{chrom_lower}_"):
                    sequence_string = f.fetch(reference=ref)
                    data[ref] = bytearray(sequence_string, 'ascii')
    return data


def run_edlib_fallback(query_seq: str, chrom: str, location: int, initial_buffer: int, max_tolerance: int,
                       error_threshold: float) -> float:
    """
    Executes an exact edit distance search via edlib, exponentially expanding the search window
    starting from initial_buffer up to max_tolerance.
    Returns the calculated error rate, exiting early if it falls below the error_threshold.
    """
    global global_sample

    target_keys = [k for k in global_sample.keys() if k == chrom or k.startswith(f"{chrom}_")]
    if not target_keys:
        return 1.0

    best_error = 1.0
    current_tolerance = initial_buffer

    prev_bounds = {k: (-1, -1) for k in target_keys}

    while True:
        expanded_any = False

        for target_key in target_keys:
            target_len = len(global_sample[target_key])

            search_start = max(0, location - current_tolerance)
            search_end = min(target_len, location + len(query_seq) + current_tolerance)

            if (search_start, search_end) == prev_bounds[target_key]:
                continue

            expanded_any = True
            prev_bounds[target_key] = (search_start, search_end)

            target_seq = global_sample[target_key][search_start:search_end].decode('ascii')

            result = edlib.align(query_seq.upper(), target_seq.upper(), mode="HW", task="path")

            if result and result['editDistance'] >= 0:
                edit_dist = result['editDistance']

                # Check that both coordinates actually exist before doing math
                if result.get('locations') and result['locations'][0][0] is not None and result['locations'][0][1] is not None:
                    loc = result['locations'][0]
                    target_match_len = loc[1] - loc[0] + 1
                    denominator = max(len(query_seq), target_match_len)
                else:
                    denominator = len(query_seq)

                error_rate = edit_dist / denominator if denominator > 0 else 1.0
                best_error = min(best_error, error_rate)

                if best_error <= error_threshold:
                    return best_error

        if not expanded_any or current_tolerance >= max_tolerance:
            break

        current_tolerance = max(1, current_tolerance * 10)
        if current_tolerance > max_tolerance:
            current_tolerance = max_tolerance

    return best_error


def score_sv(records: List[VariantRecord], buffer: Union[int, float], location_tolerance: Union[int, float],
             error_threshold=0.1):
    """
    For a single SV record, simulate the sequence and search for it in the sample genome

    :param records: List of VariantRecord objects describing the operations comprising this SV
    :param buffer: The size in bps of context to match before and after the SV
    :param location_tolerance: Allowance in bps for the SV to be matched at a different location
    :param error_threshold: The proportion of bps that can mismatch between the simulated sequence and the sample
    :return:
    """
    global global_aligners
    try:
        svid = records[0].info['SVID']
        sv_type = records[0].info['SVTYPE']

        sequences = simulate_subsequences(records, buffer)
        match_scores = []

        mappy_passed_all = True

        for query in sequences:
            sequence = query['sequence']
            chrom = query['chrom']

            best_score = 1.0

            if chrom in global_aligners:
                aligner = global_aligners[chrom]
                all_alignments = list(aligner.map(sequence))
                alignments = [a for a in all_alignments if abs(a.r_st - query['location']) <= location_tolerance]

                if len(alignments) > 0:
                    best_score = min([check_match(a, sequence) for a in alignments])

            # Run edlib fallback if mappy missed or failed the error threshold
            if len(sequence) < MIN_EDLIB_QUERY and best_score > error_threshold:
                mappy_passed_all = False
                edlib_score = run_edlib_fallback(
                    sequence,
                    chrom,
                    query['location'],
                    int(buffer),
                    int(location_tolerance),
                    error_threshold
                )
                best_score = min(best_score, edlib_score)

            match_scores.append(best_score)

        coords = records[0].pos, records[0].stop

        result = {}

        if len(sequences) > 0 and all([score <= error_threshold for score in match_scores]):
            is_correct = True
            hit_miss = 'hit'
        else:
            hit_miss = 'miss'
            is_correct = False

        rescued_by_edlib = is_correct and not mappy_passed_all

        return {
            'svid': svid,
            'sv_type': sv_type,
            'is_correct': is_correct,
            'rescued_by_edlib': rescued_by_edlib,
            'coords': coords,
            'match_scores': match_scores,
            'sequences': sequences,
        }
    except Exception as e:
        error_msg = f"Worker failed for SVID: {records[0].info.get('SVID', 'Unknown')}. Error: {e}\n{traceback.format_exc()}"
        logger.error(error_msg)
        raise e


class AlignScorer(object):
    def __init__(self, callset_vcf, buffer, gap_file):
        """
        :param callset_vcf: File containing SV callset
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
                            target_chrom and self.exclude_list[target_chrom].overlap(rec.info['TARGET'],
                                                                                     rec.info['TARGET'] + 1):
                        # record is in the excluded regions
                        allowed = False
                        break
                if allowed:
                    filtered_variants.append((svid, records))

            self.variants = {key: value for (key, value) in filtered_variants}

        self.buffer = buffer

    def load_exclude_list(self, gap_file):
        exclude_list = defaultdict(IntervalTree)
        with open(gap_file, 'r') as f:
            for line in f:
                row = line.strip().split()

                chrom = row[1]
                start, stop = int(row[2]), int(row[3])
                region_type = row[7]

                exclude_list[chrom][start:stop] = region_type

        return exclude_list

    def read_vcf(self, vcf_path: str):
        """
        Read callset .vcf and group SVs by SVID info field
        :param vcf_path: path to .vcf file
        """
        grouped_variants = defaultdict(list)
        vcf_in = pysam.VariantFile(vcf_path)

        for rec in vcf_in.fetch():
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
        edlib_rescues = 0
        precision = {}

        with ThreadPoolExecutor(max_workers=n_threads) as executor:
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

            pbar = tqdm(as_completed(futures), total=len(self.variants), desc=f'Scoring SVs', smoothing=0)

            for future in pbar:
                try:
                    result = future.result()

                    sv_type = result['sv_type']
                    svid = result['svid']
                    sequences = result['sequences']
                    coords = result['coords']
                    hit_miss = 'hit' if result['is_correct'] else 'miss'
                    match_scores = result['match_scores']
                    logger.debug(
                        f'{svid}\t{sv_type}\t{sequences[0]["chrom"]}:{coords[0]}-{coords[1]}\t{hit_miss}\t{match_scores}')

                    total_calls[sv_type] += 1
                    overall_count += 1

                    if result['is_correct']:
                        correct_calls[sv_type] += 1
                        overall_correct += 1
                        if result.get('rescued_by_edlib'):
                            edlib_rescues += 1
                    else:
                        correct_calls[sv_type] += 0

                    # Update display
                    pbar.set_description(
                        f'Scoring SVs. Current precision {overall_correct / overall_count:.2f} ({overall_correct} / {overall_count})')
                    precision[sv_type] = correct_calls[sv_type] / total_calls[sv_type]

                except Exception as exc:
                    logger.error(f'SV processing generated an exception: {exc}')

        logger.info(f"Total SVs rescued by edlib fallback: {edlib_rescues}")

        precision['ALL'] = sum(correct_calls.values()) / sum(total_calls.values())
        correct_calls['ALL'] = sum(correct_calls.values())
        total_calls['ALL'] = sum(total_calls.values())

        return precision, correct_calls, total_calls


def export_igv_session(calls, bam, classified, timestamp, igv_prefix):
    # For now, hard coded to hg19 as reference
    import xml.etree.ElementTree as ET

    output_filename = f"./igv_sessions/session_{timestamp}.xml"

    # Build XML
    session = ET.Element("Session", genome="hg19", version="8")
    resources = ET.SubElement(session, "Resources")

    if calls:
        ET.SubElement(resources, "Resource", path=igv_prefix + calls, type="vcf")
    if classified:
        ET.SubElement(resources, "Resource", path=igv_prefix + classified, type="vcf")
    if bam:
        # Standard BAM Resource
        bam_path = igv_prefix + bam
        ET.SubElement(resources, "Resource", path=bam_path, type="bam")

        bam_panel = ET.SubElement(session, "Panel", name="Alignments", height="600")

        align_track = ET.SubElement(bam_panel, "Track",
                                    clazz="org.broad.igv.sam.AlignmentTrack",
                                    displayMode="EXPANDED",
                                    id=bam_path,
                                    name=os.path.basename(bam_path),
                                    visible="true")

        ET.SubElement(align_track, "RenderOptions",
                      colorOption="READ_STRAND",
                      duplicatesOption="FILTER",
                      groupByOption="LINKED",
                      hideSmallIndels="true",
                      linkByTag="READNAME",
                      linkedReads="true",
                      smallIndelThreshold="2")

    os.makedirs(os.path.dirname(output_filename), exist_ok=True)

    tree = ET.ElementTree(session)
    tree.write(output_filename, encoding="utf-8", xml_declaration=True)

    print(f"File saved to {output_filename}")


def main():
    timestamp = datetime.datetime.now().strftime('%Y-%m-%d_%H%M%S')

    os.makedirs('./logs', exist_ok=True)
    log_filename = os.path.join("./logs", f'reconstruction_{timestamp}.log')

    file_handler = logging.FileHandler(log_filename, mode='w')
    file_handler.setLevel(logging.DEBUG)

    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)

    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s %(levelname)-8s %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        handlers=[
            file_handler,
            console_handler
        ]
    )

    parser = argparse.ArgumentParser(description='Score VCF SV calls against a reference and sample genome')
    parser.add_argument('--reference', help='Reference genome .fa file', dest='reference')
    parser.add_argument('--sample', help='Sample genome .fa file', dest='sample')
    parser.add_argument('--calls', help='VCF containing called SVs', dest='calls')
    parser.add_argument('--location_tolerance', help='BP tolerance for matching SV location', type=int,
                        default=10000000)
    parser.add_argument('--buffer', help='Subsequence context buffer', type=int, dest='buffer', default=500)
    parser.add_argument('--gap_file',
                        help='Tab-delimited file containing regions to omit (e.g., centromere and telomere)',
                        default=None)
    parser.add_argument('--config', help='Groovi call config used to infer other params', dest='config')
    parser.add_argument('--bam', help='BAM file for generating IGV config', dest='bam')
    parser.add_argument('--classified', help='VCF file of groovi-style classified breakpoints for IGV config',
                        dest='classified')
    parser.add_argument('--igv_prefix', help='Prefix for igv session paths', default='', dest='igv_prefix')
    parser.add_argument('--chrom_cache',
                        help='Directory to save/load per-chromosome MMI indices. If omitted, uses a temporary directory.',
                        default=None)
    args = parser.parse_args()

    logger.info(f'Config: {vars(args)}')

    if args.config:
        logger.info(f'Inferring params from {args.config}')
        args = update_args_from_config(args)

    global global_ref
    global global_sample
    global global_aligners

    logger.info('Initializing scorer and loading callset')
    scorer = AlignScorer(args.calls, args.buffer, args.gap_file)

    logger.info('Finding relevant chromosomes')
    chroms = set()
    for records in scorer.variants.values():
        for record in records:
            chroms.add(record.chrom)
            if 'TARGET_CHROM' in record.info:
                chroms.add(record.info['TARGET_CHROM'])
    logger.info(f'Found {len(chroms)} referenced chromosomes in callset')

    logger.info("Loading reference and sample bytearrays")
    global_ref = load_fasta_to_bytes(args.reference, chroms)
    global_sample = load_fasta_to_bytes(args.sample, chroms)
    logger.info(f"Loaded {len(global_ref)} reference chromosomes and {len(global_sample)} sample chromosomes.")

    align_params = {
        'preset': 'map-hifi',
        'k': 15,
        'w': 5,
        'best_n': 100,
        'min_cnt': 1,
        'min_dp_score': 10,
        'min_chain_score': 1,
    }

    if args.chrom_cache:
        cache_dir = args.chrom_cache
        logger.info(f"Using persistent cache directory: {cache_dir}")
    else:
        cache_dir = os.path.join(tempfile.gettempdir(), 'mappy_chrom_cache')
        logger.info(f"Using temporary cache directory: {cache_dir}")

    os.makedirs(cache_dir, exist_ok=True)

    logger.info("Building/loading per-chromosome aligners...")
    global_aligners = {}

    for chrom in chroms:
        aligner = get_chrom_aligner(args.sample, chrom, cache_dir, align_params, threads=16)
        if aligner:
            global_aligners[chrom] = aligner
        else:
            logger.warning(f"No sequences found for {chrom} in sample FASTA.")

    logger.info("Aligners ready")

    precision, correct_calls, total_calls = scorer.score_all(location_tolerance=args.location_tolerance)

    import pandas as pd

    df = pd.DataFrame({
        'correct_calls': correct_calls,
        'total_calls': total_calls,
        'precision': precision,
    })
    df.sort_values('total_calls', ascending=False, inplace=True)

    logger.info(f'Score table:\n{df}')

    export_igv_session(args.calls, args.bam, args.classified, timestamp, args.igv_prefix)


if __name__ == '__main__':
    main()