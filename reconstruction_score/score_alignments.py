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
from typing import Dict, List, Union

import pandas as pd
import pysam
from intervaltree import IntervalTree
from pysam import VariantRecord
from tqdm import tqdm

from sv_scoring_utils import (
    get_chrom_aligner, update_args_from_config, get_start_stop,
    reverse_complement, check_match, load_fasta_to_bytes, run_edlib_fallback,
    edlib_to_cigartuples, validate_junctions_from_cigar
)

logger = logging.getLogger(__name__)

MIN_EDLIB_QUERY = 5000
JUNCTION_VALIDATION_WINDOW = 150

global_ref = None
global_sample = None
global_aligners = None


def simulate_subsequences(records: List[VariantRecord], buffer: int) -> List[Dict]:
    global global_ref

    sv_type = records[0].info['SVTYPE']
    svid = records[0].info['SVID']
    chrom = records[0].chrom

    target_records = []
    in_place_records = []

    for rec in records:
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
    junction_mask = [False] * len(new_sequence)

    if not merged_target_sequence:
        new_target_sequence = list(global_ref[chrom][target_offset:target_sequence_end].decode('ascii'))
        target_changed_mask = [False] * len(new_target_sequence)
        target_junction_mask = [False] * len(new_target_sequence)

    delete_placeholder = ''
    queries = []

    def mark_junction(mask, idx):
        if 0 <= idx < len(mask):
            mask[idx] = True

    for rec in records:
        start, stop = get_start_stop(rec)
        target = rec.info.get('TARGET', stop + 1) - target_offset
        start -= offset
        stop -= offset

        if rec.info.get('TARGET_CHROM', chrom) != chrom:
            logger.warning(
                f'Skipping record: {rec.id}-{sv_type} because interchromosome target. Interchromosome checks not implemented yet.')
            return []

        if rec.info['OP_TYPE'] == 'CUT' or rec.info['SVTYPE'] == 'DEL':
            new_sequence[start:stop] = [delete_placeholder] * (stop - start)
            changed_mask[start:stop] = [True] * (stop - start)
            mark_junction(junction_mask, start)
            mark_junction(junction_mask, stop - 1)
        elif rec.info['OP_TYPE'] == 'INV' or rec.info['SVTYPE'] == 'INV':
            new_sequence[start:stop] = reverse_complement(orig_sequence[start:stop])
            changed_mask[start:stop] = [True] * (stop - start)
            mark_junction(junction_mask, start)
            mark_junction(junction_mask, stop - 1)
        elif rec.info['OP_TYPE'] == 'COPY-PASTE' or rec.info['SVTYPE'] in ['DUP', 'dDUP']:
            clip = orig_sequence[start:stop]
            if merged_target_sequence:
                new_sequence = new_sequence[:target] + clip + new_sequence[target:]
                changed_mask = changed_mask[:target] + [True] * len(clip) + changed_mask[target:]
                junction_mask = junction_mask[:target] + [False] * len(clip) + junction_mask[target:]
                mark_junction(junction_mask, target)
                mark_junction(junction_mask, target + len(clip) - 1)
            else:
                new_target_sequence = new_target_sequence[:target] + clip + new_target_sequence[target:]
                target_changed_mask = target_changed_mask[:target] + [True] * len(clip) + target_changed_mask[target:]
                target_junction_mask = target_junction_mask[:target] + [False] * len(clip) + target_junction_mask[
                    target:]
                mark_junction(target_junction_mask, target)
                mark_junction(target_junction_mask, target + len(clip) - 1)
        elif rec.info['OP_TYPE'] in ['CUT-PASTE'] or rec.info['SVTYPE'] == 'nrTRA':
            clip = orig_sequence[start:stop]
            new_sequence[start:stop] = [delete_placeholder] * (stop - start)
            changed_mask[start:stop] = [True] * len(clip)
            mark_junction(junction_mask, start)
            mark_junction(junction_mask, stop - 1)

            if merged_target_sequence:
                new_sequence = new_sequence[:target] + clip + new_sequence[target:]
                changed_mask = changed_mask[:target] + [True] * len(clip) + changed_mask[target:]
                junction_mask = junction_mask[:target] + [False] * len(clip) + junction_mask[target:]
                mark_junction(junction_mask, target)
                mark_junction(junction_mask, target + len(clip) - 1)
            else:
                new_target_sequence = new_target_sequence[:target] + clip + new_target_sequence[target:]
                target_changed_mask = target_changed_mask[:target] + [True] * len(clip) + target_changed_mask[target:]
                target_junction_mask = target_junction_mask[:target] + [False] * len(clip) + target_junction_mask[
                    target:]
                mark_junction(target_junction_mask, target)
                mark_junction(target_junction_mask, target + len(clip) - 1)
        elif rec.info['OP_TYPE'] == 'COPYinv-PASTE' or rec.info['SVTYPE'] in ['INV_dDUP']:
            clip = reverse_complement(orig_sequence[start:stop])

            if merged_target_sequence:
                new_sequence = new_sequence[:target] + clip + new_sequence[target:]
                changed_mask = changed_mask[:target] + [True] * len(clip) + changed_mask[target:]
                junction_mask = junction_mask[:target] + [False] * len(clip) + junction_mask[target:]
                mark_junction(junction_mask, target)
                mark_junction(junction_mask, target + len(clip) - 1)
            else:
                new_target_sequence = new_target_sequence[:target] + clip + new_target_sequence[target:]
                target_changed_mask = target_changed_mask[:target] + [True] * len(clip) + target_changed_mask[target:]
                target_junction_mask = target_junction_mask[:target] + [False] * len(clip) + target_junction_mask[
                    target:]
                mark_junction(target_junction_mask, target)
                mark_junction(target_junction_mask, target + len(clip) - 1)
        elif rec.info['OP_TYPE'] == 'CUTinv-PASTE' or rec.info['SVTYPE'] in ['INV_nrTRA']:
            clip = reverse_complement(orig_sequence[start:stop])
            new_sequence[start:stop] = [delete_placeholder] * (stop - start)
            changed_mask[start:stop] = [True] * len(clip)
            mark_junction(junction_mask, start)
            mark_junction(junction_mask, stop - 1)

            if merged_target_sequence:
                new_sequence = new_sequence[:target] + clip + new_sequence[target:]
                changed_mask = changed_mask[:target] + [True] * len(clip) + changed_mask[target:]
                junction_mask = junction_mask[:target] + [False] * len(clip) + junction_mask[target:]
                mark_junction(junction_mask, target)
                mark_junction(junction_mask, target + len(clip) - 1)
            else:
                new_target_sequence = new_target_sequence[:target] + clip + new_target_sequence[target:]
                target_changed_mask = target_changed_mask[:target] + [True] * len(clip) + target_changed_mask[target:]
                target_junction_mask = target_junction_mask[:target] + [False] * len(clip) + target_junction_mask[
                    target:]
                mark_junction(target_junction_mask, target)
                mark_junction(target_junction_mask, target + len(clip) - 1)
        else:
            logger.warning(f'Unknown OP_TYPE: {rec.info["OP_TYPE"]}')

    if len(changed_mask) == 0:
        logger.warning(f'No changed intervals found for {svid}-{sv_type}')
        return []

    MERGE_TOLERANCE = 5

    queries.extend(
        get_changed_subsequences(new_sequence, changed_mask, junction_mask, MERGE_TOLERANCE, offset, buffer, svid,
                                 chrom, sv_type))
    if not merged_target_sequence:
        queries.extend(
            get_changed_subsequences(new_target_sequence, target_changed_mask, target_junction_mask, MERGE_TOLERANCE,
                                     target_offset, buffer,
                                     svid, chrom, sv_type))

    return queries


def get_changed_subsequences(new_sequence, changed_mask, junction_mask, tolerance, offset, buffer, svid, chrom,
                             sv_type):
    changed_intervals = []
    queries = []

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

        seq_slice = new_sequence[adjusted_start:adjusted_stop]
        junc_slice = junction_mask[adjusted_start:adjusted_stop]

        sequence = ''
        junctions = []
        current_str_idx = 0

        for seq_char, is_junc in zip(seq_slice, junc_slice):
            if is_junc:
                junctions.append(current_str_idx)
            sequence += seq_char
            current_str_idx += len(seq_char)

        junctions = sorted(list(set(junctions)))

        query = {
            'chrom': chrom,
            'svtype': sv_type,
            'svid': svid,
            'sequence': sequence,
            'location': adjusted_start + offset,
            'length': len(sequence),
            'junctions': junctions
        }
        queries.append(query)

    return queries


def score_sv(records: List[VariantRecord], buffer: Union[int, float], location_tolerance: Union[int, float],
             error_threshold=0.1):
    global global_aligners
    global global_sample
    try:
        svid = records[0].info['SVID']
        sv_type = records[0].info['SVTYPE']

        sequences = simulate_subsequences(records, buffer)
        match_scores = []
        mappy_passed_all = True

        sv_junction_rejected = False
        best_rejected_err = None

        for query in sequences:
            sequence = query['sequence']
            chrom = query['chrom']

            best_score = 1.0
            mappy_matched = False

            seq_had_junction_rejection = False
            seq_best_rejected_err = 1.0

            if chrom in global_aligners:
                for aligner in global_aligners[chrom]:
                    all_alignments = list(aligner.map(sequence))
                    alignments = [a for a in all_alignments if abs(a.r_st - query['location']) <= location_tolerance]

                    for a in alignments:
                        err = check_match(a, sequence)
                        if err <= error_threshold:
                            if validate_junctions_from_cigar(a.cigar, query.get('junctions', []),
                                                             window=JUNCTION_VALIDATION_WINDOW,
                                                             q_st=a.q_st,
                                                             error_threshold=error_threshold):
                                best_score = min(best_score, err)
                                mappy_matched = True
                            else:
                                seq_had_junction_rejection = True
                                seq_best_rejected_err = min(seq_best_rejected_err, err)

            if len(sequence) < MIN_EDLIB_QUERY and not mappy_matched:
                mappy_passed_all = False
                for samp_bytes in global_sample:
                    edlib_res = run_edlib_fallback(
                        sequence,
                        chrom,
                        query['location'],
                        int(buffer),
                        int(location_tolerance),
                        error_threshold,
                        samp_bytes
                    )
                    if edlib_res and edlib_res['error'] <= error_threshold:
                        cigartuples = edlib_to_cigartuples(edlib_res['cigar'])
                        if validate_junctions_from_cigar(cigartuples, query.get('junctions', []),
                                                         window=JUNCTION_VALIDATION_WINDOW,
                                                         error_threshold=error_threshold):
                            best_score = min(best_score, edlib_res['error'])
                        else:
                            seq_had_junction_rejection = True
                            seq_best_rejected_err = min(seq_best_rejected_err, edlib_res['error'])

            match_scores.append(best_score)

            if best_score > error_threshold and seq_had_junction_rejection:
                sv_junction_rejected = True
                if best_rejected_err is None:
                    best_rejected_err = seq_best_rejected_err
                else:
                    best_rejected_err = min(best_rejected_err, seq_best_rejected_err)

        coords = records[0].pos, records[0].stop

        if len(sequences) > 0 and all([score <= error_threshold for score in match_scores]):
            is_correct = True
            hit_miss = 'hit'
        else:
            hit_miss = 'miss'
            is_correct = False

            if sv_junction_rejected:
                logger.debug(
                    f'Rejected {svid} {sv_type}: Alignments passed overall error (best: {best_rejected_err:.4f}) but failed junction validation.')

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
        error_msg = f'Worker failed for SVID: {records[0].info.get("SVID", "Unknown")}. Error: {e}\n{traceback.format_exc()}'
        logger.error(error_msg)
        raise e


class AlignScorer(object):
    def __init__(self, callset_vcf, buffer, gap_file):
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
        grouped_variants = defaultdict(list)
        vcf_in = pysam.VariantFile(vcf_path)
        for rec in vcf_in.fetch():
            svid = rec.info.get('SVID')
            if svid:
                grouped_variants[svid].append(rec)
        self.variants = grouped_variants
        self.vcf_header = vcf_in.header

    def score_all(self, location_tolerance=float('inf'), error_threshold=0.1, n_threads=40):
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

                    pbar.set_description(
                        f'Scoring SVs. Current precision {overall_correct / overall_count:.2f} ({overall_correct} / {overall_count})')
                    precision[sv_type] = correct_calls[sv_type] / total_calls[sv_type]

                except Exception as exc:
                    logger.error(f'SV processing generated an exception: {exc}')

        logger.info(f'Total SVs rescued by edlib fallback: {edlib_rescues}')

        precision['ALL'] = sum(correct_calls.values()) / sum(total_calls.values())
        correct_calls['ALL'] = sum(correct_calls.values())
        total_calls['ALL'] = sum(total_calls.values())

        return precision, correct_calls, total_calls


def export_igv_session(calls, bam, classified, timestamp, igv_prefix):
    import xml.etree.ElementTree as ET
    output_filename = f'./igv_sessions/session_{timestamp}.xml'

    session = ET.Element('Session', genome='hg19', version='8')
    resources = ET.SubElement(session, 'Resources')

    if calls:
        ET.SubElement(resources, 'Resource', path=igv_prefix + calls, type='vcf')
    if classified:
        ET.SubElement(resources, 'Resource', path=igv_prefix + classified, type='vcf')
    if bam:
        bam_path = igv_prefix + bam
        ET.SubElement(resources, 'Resource', path=bam_path, type='bam')
        bam_panel = ET.SubElement(session, 'Panel', name='Alignments', height='600')
        align_track = ET.SubElement(bam_panel, 'Track', clazz='org.broad.igv.sam.AlignmentTrack',
                                    displayMode='EXPANDED', id=bam_path, name=os.path.basename(bam_path),
                                    visible='true')
        ET.SubElement(align_track, 'RenderOptions', colorOption='READ_STRAND', duplicatesOption='FILTER',
                      groupByOption='LINKED', hideSmallIndels='true', linkByTag='READNAME',
                      linkedReads='true', smallIndelThreshold='2')

    os.makedirs(os.path.dirname(output_filename), exist_ok=True)
    tree = ET.ElementTree(session)
    tree.write(output_filename, encoding='utf-8', xml_declaration=True)
    print(f'File saved to {output_filename}')


def main():
    timestamp = datetime.datetime.now().strftime('%Y-%m-%d_%H%M%S')

    os.makedirs('./logs', exist_ok=True)
    log_filename = os.path.join('./logs', f'reconstruction_{timestamp}.log')

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

    parser = argparse.ArgumentParser(description='Score VCF SV calls against a reference and sample genome')
    parser.add_argument('--reference', help='Reference genome .fa file', dest='reference')
    parser.add_argument('--sample', help='Sample genome .fa file(s)', dest='sample', nargs='+')
    parser.add_argument('--calls', help='VCF containing called SVs', dest='calls')
    parser.add_argument('--location_tolerance', help='BP tolerance for matching SV location', type=int,
                        default=10000000)
    parser.add_argument('--buffer', help='Subsequence context buffer', type=int, dest='buffer', default=500)
    parser.add_argument('--gap_file', help='Tab-delimited file containing regions to omit (e.g., centromere and telomere)',
                        default=None)
    parser.add_argument('--config', help='Groovi call config used to infer other params', dest='config')
    parser.add_argument('--bam', help='BAM file for generating IGV config', dest='bam')
    parser.add_argument('--classified', help='VCF file of groovi-style classified breakpoints for IGV config',
                        dest='classified')
    parser.add_argument('--igv_prefix', help='Prefix for igv session paths', default='', dest='igv_prefix')
    parser.add_argument('--chrom_cache', help='Directory to save/load per-chromosome MMI indices.', default=None)
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

    logger.info('Loading reference and sample bytearrays')
    global_ref = load_fasta_to_bytes(args.reference, chroms)
    global_sample = [load_fasta_to_bytes(samp, chroms) for samp in args.sample]
    logger.info(f'Loaded {len(global_ref)} reference chromosomes and initialized {len(global_sample)} sample files.')

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
        logger.info(f'Using persistent cache directory: {cache_dir}')
    else:
        cache_dir = os.path.join(tempfile.gettempdir(), 'mappy_chrom_cache')
        logger.info(f'Using temporary cache directory: {cache_dir}')

    os.makedirs(cache_dir, exist_ok=True)

    logger.info('Building/loading per-chromosome aligners...')
    global_aligners = defaultdict(list)

    for samp in args.sample:
        samp_path_hash = hashlib.md5(os.path.abspath(samp).encode('utf-8')).hexdigest()[:8]
        samp_cache_dir = os.path.join(cache_dir, f'{os.path.basename(samp)}_{samp_path_hash}')
        os.makedirs(samp_cache_dir, exist_ok=True)

        for chrom in chroms:
            aligner = get_chrom_aligner(samp, chrom, samp_cache_dir, align_params, threads=32)
            if aligner:
                global_aligners[chrom].append(aligner)
            else:
                logger.warning(f'No sequences found for {chrom} in sample FASTA {samp}.')

    logger.info('Aligners ready')

    precision, correct_calls, total_calls = scorer.score_all(location_tolerance=args.location_tolerance)

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