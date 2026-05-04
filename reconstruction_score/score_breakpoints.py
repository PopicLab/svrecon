import argparse
import datetime
import logging
import os
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed

import mappy
import pysam
import pandas as pd
from tqdm import tqdm
import sys

from score_alignments import check_match, get_cached_aligner

logger = logging.getLogger(__name__)

global_aligner = None
global_ref = None


def reverse_complement(seq: str) -> str:
    """Returns the reverse complement of a given DNA sequence string."""
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
                  'a': 't', 'c': 'g', 'g': 'c', 't': 'a',
                  'N': 'N', 'n': 'n'}
    return ''.join(complement.get(base, base) for base in reversed(seq))


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
    if orientation == 'LR':
        seq = left_start + right_end
    elif orientation == 'RL':
        seq = reverse_complement(right_start) + reverse_complement(left_end)
    elif orientation == 'LL':
        seq = left_start + reverse_complement(left_end)
    elif orientation == 'RR':
        seq = reverse_complement(right_start) + right_end

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


def score_breakpoint(rec: pysam.VariantRecord, buffer: int, location_tolerance: float, error_threshold: float):
    """Simulates the junction sequence and maps it to the assembly to verify the breakpoint."""
    junction = extract_junction_sequence(rec, buffer)

    sequence = junction['sequence']
    chrom1 = junction['chrom1']
    chrom2 = junction['chrom2']
    start = junction['start']
    stop = junction['stop']
    orientation = junction['orientation']
    source = junction['source']
    svid = junction['svid']

    # Map the simulated junction sequence to the sample assembly
    all_alignments = list(global_aligner.map(sequence))

    # Filter alignments to ensure they map to one of the expected chromosomes
    alignments = [a for a in all_alignments if a.ctg.startswith(chrom1) or a.ctg.startswith(chrom2)]

    # Filter alignments by distance from the expected breakpoint coordinate if a tolerance is provided
    if location_tolerance != float('inf'):
        alignments = [a for a in alignments if
                      abs(a.r_st - start) <= location_tolerance or abs(a.r_st - stop) <= location_tolerance]

    # Calculate the match error rate and determine if it falls within the acceptable threshold
    if len(alignments) > 0:
        match_scores = [check_match(a, sequence) for a in alignments]
        min_error = min(match_scores)
        is_correct = min_error <= error_threshold
    else:
        is_correct = False
        match_scores = [1.0]

    return {
        'svid': svid,
        'orientation': orientation,
        'source': source,
        'is_correct': is_correct,
        'match_scores': match_scores,
        'coords': (start, stop),
        'chrom': chrom1
    }


def main():
    """Orchestrates the loading of data, parallel execution of scoring, and aggregation of results."""
    # Setup logging configuration
    timestamp = datetime.datetime.now().strftime('%Y-%m-%d_%H%M%S')
    os.makedirs("./logs", exist_ok=True)
    log_filename = os.path.join("./logs", f'score_breakpoints_{timestamp}.log')

    # Create handlers and set their individual log levels
    file_handler = logging.FileHandler(log_filename, mode='w')
    file_handler.setLevel(logging.DEBUG)  # File gets DEBUG and above

    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)  # Console gets INFO and above (skips DEBUG)

    # Pass the configured handlers to basicConfig
    logging.basicConfig(
        level=logging.DEBUG,  # The root logger must be set to the lowest level you want to capture
        format='%(asctime)s %(levelname)-8s %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        handlers=[
            file_handler,
            console_handler
        ]
    )

    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description='Score VCF breakpoints against a sample assembly via simulated junctions')
    parser.add_argument('--reference', required=True, help='Reference genome .fa file')
    parser.add_argument('--sample', required=True, help='Sample genome/assembly .fa file')
    parser.add_argument('--calls', required=True, help='VCF containing called SV breakpoints')
    parser.add_argument('--buffer', help='Subsequence context buffer length', type=int, default=50)
    parser.add_argument('--location_tolerance', help='Max distance for alignment coordinate match', type=float,
                        default=float('inf'))
    parser.add_argument('--error_threshold', help='Max sequence error rate', type=float, default=0.1)
    parser.add_argument('--threads', help='Number of threads', type=int, default=16)
    parser.add_argument('--min_support', help='Minimum support to score a breakpoint', type=float, default=1)
    args = parser.parse_args()

    logger.info(f'Config: {vars(args)}')

    global global_ref
    global global_aligner

    # Parse the input VCF file and load records into memory
    logger.info('Parsing VCF callset...')
    vcf_in = pysam.VariantFile(args.calls)
    records = list(vcf_in.fetch())
    logger.info(f'Loaded {len(records)} breakpoint records.')

    # Identify all unique chromosomes present in the breakpoint records
    logger.info('Finding relevant chromosomes...')
    chroms = set()
    for rec in records:
        chroms.add(rec.chrom)
        chroms.add(rec.info['CHR2'])
    logger.info(f'Found {len(chroms)} referenced chromosomes in callset.')

    # Load only the required chromosome sequences from the reference fasta into a bytearray dictionary
    logger.info("Loading reference into memory...")
    ref_data = {}
    with pysam.FastaFile(args.reference) as f:
        for chrom in chroms:
            sequence_string = f.fetch(reference=chrom)
            ref_data[chrom] = bytearray(sequence_string, 'ascii')
    global_ref = ref_data
    logger.info(f"Reference pre-loaded with {len(global_ref)} chromosomes.")

    # Initialize the mappy aligner with the sample assembly fasta
    logger.info("Loading index into aligner...")

    global_aligner = get_cached_aligner(fasta_path=args.sample, threads=args.threads, preset='map-pb')

    logger.info("Aligner ready.")

    total_calls = Counter()
    correct_calls = Counter()
    overall_count = 0
    overall_correct = 0

    filtered_records = [rec for rec in records if rec.info['SUPPORT'] >= args.min_support]

    # Execute the breakpoint scoring in parallel across the specified number of threads
    logger.info("Scoring Breakpoints...")
    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        futures = {
            executor.submit(
                score_breakpoint,
                rec,
                args.buffer,
                args.location_tolerance,
                args.error_threshold,
            )
            for rec in filtered_records
        }

        pbar = tqdm(as_completed(futures), total=len(filtered_records), desc='Scoring Breakpoints')

        # Aggregate the results as each thread completes
        for future in pbar:
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

            if overall_count > 0:
                pbar.set_description(
                    f'Scoring. Precision {overall_correct / overall_count:.2f} ({overall_correct}/{overall_count})')

            logger.debug(
                f"{svid}\t{source}\t{orientation}\t{result['chrom']}:{result['coords'][0]}-{result['coords'][1]}\t{hit_miss}\t{result['match_scores']}")

    # Compile the final statistics broken down by Source and Orientation
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

    # Format and display the output dataframe
    df = pd.DataFrame(results_data)
    df.sort_values(by=['Source', 'Orientation'], ascending=[True, True], inplace=True)
    df.reset_index(drop=True, inplace=True)

    logger.info("\n--- Final Results ---")
    logger.info(f'Score table:\n{df.to_string(index=False)}')


if __name__ == '__main__':
    main()