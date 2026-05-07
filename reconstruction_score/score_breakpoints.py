import argparse
import datetime
import logging
import os
import sys
import tempfile
import traceback
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd
import pysam
from tqdm import tqdm

from sv_scoring_utils import (
    get_chrom_aligner, update_args_from_config, reverse_complement,
    check_match, load_fasta_to_bytes, run_edlib_fallback
)

logger = logging.getLogger(__name__)

global_ref = None
global_sample = None
global_aligners = None


def get_sequence_slice(chrom: str, start: int, stop: int) -> str:
    """Extracts a specific subsequence from the pre-loaded global reference bytearray."""
    return global_ref[chrom][start:stop].decode('ascii')

def extract_junction_sequence(rec: pysam.VariantRecord, buffer: int) -> dict:
    """Extracts the expected sample sequence across the breakpoint junction based on ORIENTATION rules."""
    chrom1 = rec.chrom
    start = rec.pos - 1

    info = rec.info
    stop = rec.stop - 1
    chrom2 = info['CHR2']

    # Pysam may parse INFO fields as tuples depending on the VCF header
    orientation = info['ORIENTATION']
    if isinstance(orientation, tuple):
        orientation = orientation[0]
    if isinstance(orientation, str) and orientation.startswith('k'):
        orientation = orientation[1:]

    source = info['SOURCE']
    if isinstance(source, tuple):
        source = source[0]
    if isinstance(source, str) and source.startswith('k'):
        source = source[1:]

    if orientation is None and source is None:
        orientation = info.get('BKP_TYPE', '')
        source = info.get('BKP_TYPE', '')

    # Extract the flanking sequences around the start and stop breakpoints
    left_start = get_sequence_slice(chrom1, start - buffer, start)
    right_start = get_sequence_slice(chrom1, start, start + buffer)

    left_end = get_sequence_slice(chrom2, stop - buffer, stop)
    right_end = get_sequence_slice(chrom2, stop, stop + buffer)

    # Build the sequence across the junction based on the orientation
    # Cast the list output of reverse_complement back to strings
    if orientation == 'LR':
        seq = left_start + right_end
    elif orientation == 'RL':
        seq = "".join(reverse_complement(right_start)) + "".join(reverse_complement(left_end))
    elif orientation == 'LL':
        seq = left_start + "".join(reverse_complement(left_end))
    elif orientation == 'RR':
        seq = "".join(reverse_complement(right_start)) + right_end
    else:
        seq = ""

    return {
        'chrom1': chrom1,
        'chrom2': chrom2,
        'start': start,
        'stop': stop,
        'orientation': orientation,
        'source': source,
        'sequence': seq,
        'svid': rec.id or f"{chrom1}_{start}_{orientation}",
        'support': info['SUPPORT']
    }


def score_breakpoint(rec: pysam.VariantRecord, buffer: int, location_tolerance: int, error_threshold: float):
    """Simulates the junction sequence and maps it to the assembly to verify the breakpoint."""
    global global_aligners
    global global_sample

    junction = extract_junction_sequence(rec, buffer)

    sequence = junction['sequence']
    chrom1 = junction['chrom1']
    chrom2 = junction['chrom2']
    start = junction['start']
    stop = junction['stop']
    orientation = junction['orientation']
    source = junction['source']
    svid = junction['svid']

    mappy_passed_all = True
    best_score = 1.0

    # Search against the aligners for BOTH chromosomes involved in the junction
    target_chroms = {chrom1, chrom2}
    alignments = []

    for tgt_chrom in target_chroms:
        if tgt_chrom in global_aligners:
            for a in global_aligners[tgt_chrom].map(sequence):
                # Ensure the mapped contig matches one of the expected chromosomes
                if a.ctg.startswith(chrom1) or a.ctg.startswith(chrom2):
                    if abs(a.r_st - start) <= location_tolerance or abs(a.r_st - stop) <= location_tolerance:
                        alignments.append(a)

    if len(alignments) > 0:
        best_score = min([check_match(a, sequence) for a in alignments])

    # Run edlib fallback if mappy missed
    if best_score > error_threshold:
        mappy_passed_all = False

        # Edlib fallback on chrom1 around the origin breakpoint
        edlib_score1 = run_edlib_fallback(
            sequence, chrom1, start, int(buffer),
            location_tolerance, error_threshold, global_sample
        )
        # Edlib fallback on chrom2 around the target breakpoint
        edlib_score2 = run_edlib_fallback(
            sequence, chrom2, stop, int(buffer),
            location_tolerance, error_threshold, global_sample
        )
        best_score = min(best_score, edlib_score1, edlib_score2)

    is_correct = best_score <= error_threshold
    rescued_by_edlib = is_correct and not mappy_passed_all

    return {
        'svid': svid,
        'orientation': orientation,
        'source': source,
        'is_correct': is_correct,
        'rescued_by_edlib': rescued_by_edlib,
        'match_scores': [best_score],
        'coords': (start, stop),
        'chrom': chrom1
    }


def main():
    """Orchestrates the loading of data, parallel execution of scoring, and aggregation of results."""
    timestamp = datetime.datetime.now().strftime('%Y-%m-%d_%H%M%S')
    os.makedirs("./logs", exist_ok=True)
    log_filename = os.path.join("./logs", f'score_breakpoints_{timestamp}.log')

    file_handler = logging.FileHandler(log_filename, mode='w')
    file_handler.setLevel(logging.DEBUG)

    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)

    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s %(levelname)-8s %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        handlers=[file_handler, console_handler]
    )

    parser = argparse.ArgumentParser(description='Score VCF breakpoints against a sample assembly via simulated junctions')
    parser.add_argument('--reference', help='Reference genome .fa file')
    parser.add_argument('--sample', help='Sample genome/assembly .fa file')
    parser.add_argument('--calls', help='VCF containing called SV breakpoints')
    parser.add_argument('--config', help='Groovi call config used to infer other params', dest='config')
    parser.add_argument('--buffer', help='Subsequence context buffer length', type=int, default=50)
    parser.add_argument('--location_tolerance', help='Max distance for alignment coordinate match', type=int, default=10000000)
    parser.add_argument('--error_threshold', help='Max sequence error rate', type=float, default=0.1)
    parser.add_argument('--threads', help='Number of threads', type=int, default=16)
    parser.add_argument('--min_support', help='Minimum support to score a breakpoint', type=float, default=1)
    parser.add_argument('--chrom_cache', help='Directory to save/load per-chromosome MMI indices.', default=None)
    args = parser.parse_args()

    logger.info(f'Config: {vars(args)}')

    if args.config:
        logger.info(f'Inferring params from {args.config}')
        args = update_args_from_config(args)

    global global_ref
    global global_sample
    global global_aligners

    logger.info('Parsing VCF callset...')
    vcf_in = pysam.VariantFile(args.calls)
    records = list(vcf_in.fetch())
    logger.info(f'Loaded {len(records)} breakpoint records.')

    logger.info('Finding relevant chromosomes...')
    chroms = set()
    for rec in records:
        chroms.add(rec.chrom)
        chroms.add(rec.info['CHR2'])
    logger.info(f'Found {len(chroms)} referenced chromosomes in callset.')

    logger.info("Loading reference and sample bytearrays")
    global_ref = load_fasta_to_bytes(args.reference, chroms)
    global_sample = load_fasta_to_bytes(args.sample, chroms)

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
    else:
        cache_dir = os.path.join(tempfile.gettempdir(), 'mappy_chrom_cache')
    os.makedirs(cache_dir, exist_ok=True)

    logger.info("Building/loading per-chromosome aligners...")
    global_aligners = {}
    for chrom in chroms:
        aligner = get_chrom_aligner(args.sample, chrom, cache_dir, align_params, threads=16)
        if aligner:
            global_aligners[chrom] = aligner

    logger.info("Aligners ready.")

    total_calls = Counter()
    correct_calls = Counter()
    overall_count = 0
    overall_correct = 0
    edlib_rescues = 0

    filtered_records = [rec for rec in records if rec.info['SUPPORT'] >= args.min_support]

    logger.info("Scoring Breakpoints...")
    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        futures = {
            executor.submit(
                score_breakpoint, rec, args.buffer, args.location_tolerance, args.error_threshold
            )
            for rec in filtered_records
        }

        pbar = tqdm(as_completed(futures), total=len(filtered_records), desc='Scoring Breakpoints', smoothing=0)

        for future in pbar:
            try:
                result = future.result()
                orientation = result['orientation']
                source = result['source']
                svid = result['svid']
                hit_miss = 'hit' if result['is_correct'] else 'miss'

                key = (source, orientation)
                total_calls[key] += 1
                overall_count += 1

                if result['is_correct']:
                    correct_calls[key] += 1
                    overall_correct += 1
                    if result.get('rescued_by_edlib'):
                        edlib_rescues += 1

                if overall_count > 0:
                    pbar.set_description(
                        f'Scoring. Precision {overall_correct / overall_count:.2f} ({overall_correct}/{overall_count})')

                logger.debug(
                    f"{svid}\t{source}\t{orientation}\t{result['chrom']}:{result['coords'][0]}-{result['coords'][1]}\t{hit_miss}\t{result['match_scores']}")
            except Exception as exc:
                logger.error(f'Breakpoint processing generated an exception: {exc}')

    logger.info(f"Total Breakpoints rescued by edlib fallback: {edlib_rescues}")

    results_data = []
    for key in total_calls:
        source, orientation = key
        tot = total_calls[key]
        cor = correct_calls[key]
        results_data.append({
            'Source': source,
            'Orientation': orientation,
            'correct_calls': cor,
            'total_calls': tot,
            'precision': cor / tot if tot > 0 else 0.0
        })

    if overall_count > 0:
        results_data.append({
            'Source': 'ALL',
            'Orientation': 'ALL',
            'correct_calls': overall_correct,
            'total_calls': overall_count,
            'precision': overall_correct / overall_count
        })

    df = pd.DataFrame(results_data)
    df.sort_values(by=['Source', 'Orientation'], ascending=[True, True], inplace=True)
    df.reset_index(drop=True, inplace=True)

    logger.info("\n--- Final Results ---")
    logger.info(f'Score table:\n{df.to_string(index=False)}')


if __name__ == '__main__':
    main()