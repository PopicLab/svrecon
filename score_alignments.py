import argparse
import datetime
import hashlib
import json
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
    edlib_to_cigartuples, validate_junctions_from_cigar,
    BamReader, run_read_edlib
)

logger = logging.getLogger(__name__)

MIN_EDLIB_QUERY = 5000
# Search-window radius for the assembly edlib fallback (a legacy path for short
# queries mappy missed). This is a COST bound on an unbounded O(n*m) edlib search,
# kept separate from --location_tolerance: the latter only filters mappy hits by
# genomic position and defaults to infinity (meaningless for unscaffolded, per-
# contig assemblies where a hit's r_st is contig-local, not a genomic coordinate).
EDLIB_FALLBACK_MAX_TOLERANCE = 10_000_000
JUNCTION_VALIDATION_WINDOW = 300

global_ref = None
global_sample = None
global_aligners = None

# Read-based evaluation config (set in main()).
global_eval_mode = 'assembly'          # 'assembly' | 'reads' | 'both'
global_bam_reader = None
global_read_error_threshold = 0.1
global_min_read_support = 1
global_max_reads_per_site = 1000  # overridden by --max_reads_per_site in main()


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
            clip = list(reverse_complement(orig_sequence[start:stop]))

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
            clip = list(reverse_complement(orig_sequence[start:stop]))
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
    # Length of the FULL resulting allele these subsequences are carved from -- a
    # read must be at least this long to contain the whole variant (read mode).
    # Sum the element lengths rather than len(new_sequence): deletions leave empty
    # '' placeholders (and inserted clips can be multi-char), so element count
    # over-estimates the true bp length for deletion-containing alleles.
    result_len = sum(len(s) for s in new_sequence)

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
            'junctions': junctions,
            'result_len': result_len,
        }
        queries.append(query)

    return queries


def score_sv(records: List[VariantRecord], buffer: Union[int, float], location_tolerance: Union[int, float],
             error_threshold=0.1):
    global global_aligners
    global global_sample
    global global_eval_mode, global_bam_reader, global_read_error_threshold
    global global_min_read_support, global_max_reads_per_site
    try:
        svid = records[0].info['SVID']
        sv_type = records[0].info['SVTYPE']

        sequences = simulate_subsequences(records, buffer)
        match_scores = []
        subseq_diagnostics = []
        subseq_sources = []
        subseq_status = []   # per-subsequence: 'pass' | 'fail' | 'inconclusive'
        segments = []        # structured per-subsequence detail for the optional JSON report
        mappy_passed_all = True

        sv_junction_rejected = False
        best_rejected_err = None

        # Reference footprint of the SV (incl. any DUP/translocation TARGET loci),
        # used as the window reads must OVERLAP to become candidates in
        # read-based evaluation (we don't require spanning; SV-carrying reads are
        # typically split-aligned).
        read_bp_start = read_bp_end = None
        if global_eval_mode in ('reads', 'both'):
            starts = [r.start for r in records]
            stops = [r.stop for r in records]
            for r in records:
                if 'TARGET' in r.info:
                    starts.append(int(r.info['TARGET']))
                    stops.append(int(r.info['TARGET']))
            read_bp_start = max(0, min(starts))
            read_bp_end = max(stops)

        for query in sequences:
            sequence = query['sequence']
            chrom = query['chrom']

            best_score = 1.0
            mappy_matched = False

            seq_had_junction_rejection = False
            seq_best_rejected_err = 1.0

            seq_has_aligner = chrom in global_aligners
            seq_saw_candidate = False
            seq_best_fail_err = None
            seq_no_reads = False
            seq_inconclusive = False  # reads overlap but none span the full resulting allele
            seq_reads_aborted = False  # spanning reads existed but all exceeded the 2x edlib bound
            seq_source = None  # which tier validated this subsequence: assembly | edlib | reads
            seq_passed = False  # a tier validated it within ITS OWN threshold (+ junctions)
            seq_junctions = []  # per-junction results from the tier whose check we report

            # Read-based evaluation
            if global_eval_mode in ('reads', 'both'):
                read_seqs = global_bam_reader.candidate_read_seqs(
                    chrom, read_bp_start, read_bp_end, int(buffer), global_max_reads_per_site)
                # A read can only confirm the variant if it is long enough to hold
                # the resulting allele
                result_len = query.get('result_len') or len(sequence)
                spanning = [r for r in read_seqs if len(r) >= result_len]
                if not read_seqs:
                    seq_no_reads = True
                elif not spanning:
                    seq_inconclusive = True
                else:
                    read_res, _ = run_read_edlib(sequence, spanning,
                                                 global_read_error_threshold, global_min_read_support)
                    if read_res:
                        seq_saw_candidate = True
                        if read_res['error'] <= global_read_error_threshold:
                            cigartuples = edlib_to_cigartuples(read_res['cigar'])
                            junc = validate_junctions_from_cigar(cigartuples, query.get('junctions', []),
                                                             window=JUNCTION_VALIDATION_WINDOW,
                                                             error_threshold=global_read_error_threshold)
                            if not seq_passed:
                                seq_junctions = junc
                            if all(j['passed'] for j in junc):
                                best_score = min(best_score, read_res['error'])
                                seq_source = 'reads'
                                seq_passed = True
                            else:
                                seq_had_junction_rejection = True
                                seq_best_rejected_err = min(seq_best_rejected_err, read_res['error'])
                                seq_best_fail_err = read_res['error'] if seq_best_fail_err is None else min(seq_best_fail_err, read_res['error'])
                        else:
                            seq_best_fail_err = read_res['error'] if seq_best_fail_err is None else min(seq_best_fail_err, read_res['error'])
                    else:
                        # spanning reads existed but all exceeded the 2x edlib bound
                        seq_reads_aborted = True

            if global_eval_mode in ('assembly', 'both') and not seq_passed and chrom in global_aligners:
                # Assembly check (Tier 2). In 'both' mode this only runs when the
                # reads above did not already validate the subsequence. Align the
                # reconstructed subsequence to the per-chromosome assembly
                # aligner(s), keep alignments within location_tolerance of the
                # expected locus, and accept a match only when the error rate
                # clears error_threshold AND the SV junctions validate.
                for aligner in global_aligners[chrom]:
                    all_alignments = list(aligner.map(sequence))
                    alignments = [a for a in all_alignments if abs(a.r_st - query['location']) <= location_tolerance]

                    for a in alignments:
                        seq_saw_candidate = True
                        err = check_match(a, sequence)
                        if err <= error_threshold:
                            junc = validate_junctions_from_cigar(a.cigar, query.get('junctions', []),
                                                             window=JUNCTION_VALIDATION_WINDOW,
                                                             q_st=a.q_st,
                                                             error_threshold=error_threshold)
                            if not seq_passed:
                                seq_junctions = junc
                            if all(j['passed'] for j in junc):
                                best_score = min(best_score, err)
                                mappy_matched = True
                                seq_source = 'assembly'
                                seq_passed = True
                            else:
                                seq_had_junction_rejection = True
                                seq_best_rejected_err = min(seq_best_rejected_err, err)
                                seq_best_fail_err = err if seq_best_fail_err is None else min(seq_best_fail_err, err)
                        else:
                            seq_best_fail_err = err if seq_best_fail_err is None else min(seq_best_fail_err, err)

            if global_eval_mode in ('assembly', 'both') and not seq_passed and len(sequence) < MIN_EDLIB_QUERY and not mappy_matched:
                # Assembly edlib fallback: for short subsequences that mappy
                # failed to align (mappy can miss very short queries), retry with
                # edlib against expanding windows of the assembly, applying the
                # same error-threshold + junction-validation acceptance criteria.
                mappy_passed_all = False
                for samp_bytes in global_sample:
                    edlib_res = run_edlib_fallback(
                        sequence,
                        chrom,
                        query['location'],
                        int(buffer),
                        EDLIB_FALLBACK_MAX_TOLERANCE,
                        error_threshold,
                        samp_bytes
                    )
                    if edlib_res:
                        seq_saw_candidate = True
                        if edlib_res['error'] <= error_threshold:
                            cigartuples = edlib_to_cigartuples(edlib_res['cigar'])
                            junc = validate_junctions_from_cigar(cigartuples, query.get('junctions', []),
                                                             window=JUNCTION_VALIDATION_WINDOW,
                                                             error_threshold=error_threshold)
                            if not seq_passed:
                                seq_junctions = junc
                            if all(j['passed'] for j in junc):
                                best_score = min(best_score, edlib_res['error'])
                                seq_source = 'edlib'
                                seq_passed = True
                            else:
                                seq_had_junction_rejection = True
                                seq_best_rejected_err = min(seq_best_rejected_err, edlib_res['error'])
                                seq_best_fail_err = edlib_res['error'] if seq_best_fail_err is None else min(seq_best_fail_err, edlib_res['error'])
                        else:
                            seq_best_fail_err = edlib_res['error'] if seq_best_fail_err is None else min(seq_best_fail_err, edlib_res['error'])

            match_scores.append(best_score)

            # Classify the subsequence from the evidence gathered above. It falls
            # into exactly one category, and the category fixes both the scoring
            # status and the diagnostic reason. Checked in priority order:
            #   CONFIRMED    -- some tier matched within threshold -> pass.
            #   CONTRADICTED -- a candidate (spanning read or assembly alignment)
            #                   was evaluated and disagreed -> fail. In 'both' mode
            #                   this means neither tier could confirm it.
            #   UNTESTABLE   -- the read tier had no read to judge with (no
            #                   overlapping / no spanning read) -> inconclusive:
            #                   absence of a read is not evidence against the call.
            #   NOT_FOUND    -- assembly tier produced nothing to compare against.
            #                   In assembly mode absence IS the verdict -> fail.
            if seq_passed:
                status, reason = 'pass', 'pass'
            elif seq_had_junction_rejection or seq_saw_candidate or seq_reads_aborted:
                status = 'fail'
                reason = 'junction_failed' if seq_had_junction_rejection else 'over_error_threshold'
            elif seq_no_reads or seq_inconclusive:
                status = 'inconclusive'
                reason = 'no_reads' if seq_no_reads else 'inconclusive'
            else:
                status = 'fail'
                reason = 'no_aligner' if not seq_has_aligner else 'no_alignment_in_window'

            if seq_best_fail_err is not None:
                err_str = f'{seq_best_fail_err:.4f}'
            elif reason == 'over_error_threshold':
                err_str = 'aborted'  # spanning reads exceeded the 2x-threshold edlib bound
            else:
                err_str = 'NA'
            subseq_diagnostics.append(f'{reason}:{err_str}')
            subseq_sources.append(seq_source)
            subseq_status.append(status)

            # Structured per-subsequence record for the optional JSON report. `error` is the
            # winning error when the subsequence passed, else the best failing error, else None
            # (nothing testable produced a number).
            if status == 'pass':
                seg_error = float(f'{best_score:.4g}')
            elif seq_best_fail_err is not None:
                seg_error = float(f'{seq_best_fail_err:.4g}')
            else:
                seg_error = None
            segments.append({
                'chrom': chrom,
                'ref_start': query['location'],
                'status': status,
                'reason': reason,
                'source': seq_source,
                'error': seg_error,
                'junctions': seq_junctions,
            })

            if best_score > error_threshold and seq_had_junction_rejection:
                sv_junction_rejected = True
                if best_rejected_err is None:
                    best_rejected_err = seq_best_rejected_err
                else:
                    best_rejected_err = min(best_rejected_err, seq_best_rejected_err)

        coords = records[0].pos, records[0].stop

        # SV-level roll-up:
        #   hit          -> every subsequence passed
        #   miss         -> at least one subsequence was contradicted (spanning
        #                   reads / assembly aligned but didn't match)
        #   inconclusive -> no contradiction, but at least one subsequence could
        #                   not be tested (reads mode: no read spans the full allele)
        if len(subseq_status) > 0 and all(s == 'pass' for s in subseq_status):
            outcome = 'hit'
        elif 'fail' in subseq_status:
            outcome = 'miss'
        elif 'inconclusive' in subseq_status:
            outcome = 'inconclusive'
        else:
            outcome = 'miss'
        is_correct = (outcome == 'hit')

        if outcome == 'miss' and sv_junction_rejected:
            logger.debug(
                f'Rejected {svid} {sv_type}: Alignments passed overall error (best: {best_rejected_err:.4f}) but failed junction validation.')

        rescued_by_edlib = is_correct and not mappy_passed_all

        # Per-SV provenance: the highest tier any passing subsequence needed.
        # 'reads' implies the SV would have been an assembly-miss without reads.
        if is_correct:
            srcs = set(s for s in subseq_sources if s)
            if 'reads' in srcs:
                validation_source = 'reads'
            elif 'edlib' in srcs:
                validation_source = 'edlib'
            else:
                validation_source = 'assembly'
        else:
            validation_source = None

        # Coarse evaluation tier for two-tier reporting:
        #   'read'     (Tier 1) -- a single actual read contained the allele. The
        #              strongest, most concrete evidence; independent of any
        #              assembly's quality.
        #   'assembly' (Tier 2) -- confirmed only against the reconstructed
        #              assembly (mappy or the edlib fallback). Supporting evidence
        #              for events too long for any single read to span.
        if validation_source == 'reads':
            tier = 'read'
        elif validation_source in ('assembly', 'edlib'):
            tier = 'assembly'
        else:
            tier = None

        return {
            'svid': svid,
            'sv_type': sv_type,
            'is_correct': is_correct,
            'outcome': outcome,   # 'hit' | 'miss' | 'inconclusive'
            'rescued_by_edlib': rescued_by_edlib,
            'validation_source': validation_source,
            'tier': tier,   # 'read' (Tier 1) | 'assembly' (Tier 2) | None
            'coords': coords,
            'match_scores': match_scores,
            'sequences': sequences,
            'subseq_diagnostics': subseq_diagnostics,
            'segments': segments,
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

    def score_all(self, location_tolerance=float('inf'), error_threshold=0.1, n_threads=40, report_path=None):
        total_calls = Counter()
        correct_calls = Counter()
        inconclusive_calls = Counter()  # reads mode: no read spans the full resulting allele
        assembly_hits = Counter()   # hits validated by the assembly (mappy or edlib)
        read_hits = Counter()       # hits validated only by reads
        source_counts = Counter()   # overall hit provenance: assembly | edlib | reads
        overall_count = 0
        overall_correct = 0
        overall_inconclusive = 0
        edlib_rescues = 0
        precision = {}

        # Optional per-SV JSON sidecar (one record per line). Additive: it does not affect the
        # log output or the score table. Written from this single consumer thread, so no lock.
        report_fh = open(report_path, 'w') if report_path else None

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
                    outcome = result.get('outcome', 'hit' if result['is_correct'] else 'miss')
                    match_scores = result['match_scores']
                    line = f'{svid}\t{sv_type}\t{sequences[0]["chrom"]}:{coords[0]}-{coords[1]}\t{outcome}\t{match_scores}'
                    if outcome == 'hit':
                        # Tier 1 (read) vs Tier 2 (assembly), plus the detailed source.
                        line += f'\t{result.get("tier", "")}\t{result.get("validation_source", "")}'
                    else:  # miss or inconclusive -> show why
                        line += f'\t{result.get("subseq_diagnostics", [])}'
                    logger.debug(line)

                    if report_fh is not None:
                        report_fh.write(json.dumps({
                            'svid': svid,
                            'svtype': sv_type,
                            'outcome': outcome,
                            'tier': result.get('tier'),
                            'segments': result.get('segments', []),
                        }) + '\n')

                    total_calls[sv_type] += 1
                    overall_count += 1

                    if outcome == 'hit':
                        correct_calls[sv_type] += 1
                        overall_correct += 1
                        if result.get('rescued_by_edlib'):
                            edlib_rescues += 1
                        src = result.get('validation_source') or 'assembly'
                        source_counts[src] += 1
                        if src == 'reads':
                            read_hits[sv_type] += 1
                        else:
                            assembly_hits[sv_type] += 1
                    elif outcome == 'inconclusive':
                        inconclusive_calls[sv_type] += 1
                        overall_inconclusive += 1
                    # else: miss -> counts toward total but not correct/inconclusive

                    # precision is over CONCLUSIVE calls only (hits + misses)
                    conclusive = overall_count - overall_inconclusive
                    pbar.set_description(
                        f'Scoring SVs. Precision {overall_correct / conclusive:.2f} '
                        f'({overall_correct}/{conclusive}); inconclusive {overall_inconclusive}')
                    denom = total_calls[sv_type] - inconclusive_calls[sv_type]
                    precision[sv_type] = correct_calls[sv_type] / denom if denom else float('nan')

                except Exception as exc:
                    logger.error(f'SV processing generated an exception: {exc}')

        if report_fh is not None:
            report_fh.close()
            logger.info(f'Wrote per-SV eval report: {report_path}')

        logger.info(f'Total SVs rescued by edlib fallback: {edlib_rescues}')
        overall_miss = overall_count - overall_correct - overall_inconclusive
        logger.info(
            f"Outcomes: hit={overall_correct} miss={overall_miss} inconclusive={overall_inconclusive} "
            f"(inconclusive = no read spans the full resulting allele; excluded from precision)")
        logger.info(
            f"Hits by validation source: assembly={source_counts['assembly']} "
            f"edlib={source_counts['edlib']} reads={source_counts['reads']}")
        logger.info(
            f"Hits by tier: read(Tier1)={source_counts['reads']} "
            f"assembly(Tier2)={source_counts['assembly'] + source_counts['edlib']} "
            f"(Tier1 = a single real read contained the allele; Tier2 = confirmed only against the assembly)")

        total_all = sum(total_calls.values())
        inc_all = sum(inconclusive_calls.values())
        conclusive_all = total_all - inc_all
        precision['ALL'] = sum(correct_calls.values()) / conclusive_all if conclusive_all else float('nan')
        correct_calls['ALL'] = sum(correct_calls.values())
        total_calls['ALL'] = total_all
        inconclusive_calls['ALL'] = inc_all
        assembly_hits['ALL'] = sum(assembly_hits.values())
        read_hits['ALL'] = sum(read_hits.values())

        return precision, correct_calls, total_calls, inconclusive_calls, assembly_hits, read_hits


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


def log_name_from_config(config):
    """Derive a stable log basename from the config's experiment directory.

    The experiment is identified by the config's directory relative to an
    `experiments/` ancestor, with path components joined by dots. For example:
        .../experiments/HG00733/groovi.yaml            -> reconstruction.HG00733
        .../experiments/hg002/latest_filtering/groovi.yaml
                                                       -> reconstruction.hg002.latest_filtering
    Falls back to the config's parent directory name when no `experiments/`
    ancestor is present, and returns None if no config is given.
    """
    if not config:
        return None

    exp_dir = os.path.dirname(os.path.abspath(config))
    parts = exp_dir.split(os.sep)
    if 'experiments' in parts:
        idx = parts.index('experiments')
        rel_parts = parts[idx + 1:]
    else:
        rel_parts = parts[-1:]

    rel_parts = [p for p in rel_parts if p]
    if not rel_parts:
        return None
    return 'reconstruction.' + '.'.join(rel_parts)


def main():
    timestamp = datetime.datetime.now().strftime('%Y-%m-%d_%H%M%S')

    parser = argparse.ArgumentParser(description='Score VCF SV calls against a reference and sample genome')
    parser.add_argument('--reference', help='Reference genome .fa file', dest='reference')
    parser.add_argument('--sample', help='Sample genome .fa file(s)', dest='sample', nargs='+')
    parser.add_argument('--calls', help='VCF containing called SVs', dest='calls')
    parser.add_argument('--location_tolerance',
                        help='BP tolerance between a mappy hit and the expected SV location; '
                             'default unbounded (recommended for per-contig/unscaffolded assemblies '
                             'whose hit coordinates are contig-local, not genomic)',
                        type=float, default=float('inf'))
    parser.add_argument('--buffer', help='Subsequence context buffer', type=int, dest='buffer', default=500)
    parser.add_argument('--gap_file', help='Tab-delimited file containing regions to omit (e.g., centromere and telomere)',
                        default=None)
    parser.add_argument('--config', help='Groovi call config used to infer other params', dest='config')
    parser.add_argument('--bam', help='BAM file for generating IGV config', dest='bam')
    parser.add_argument('--classified', help='VCF file of groovi-style classified breakpoints for IGV config',
                        dest='classified')
    parser.add_argument('--igv_prefix', help='Prefix for igv session paths', default='', dest='igv_prefix')
    parser.add_argument('--chrom_cache', help='Directory to save/load per-chromosome MMI indices.', default=None)
    parser.add_argument('--eval_mode', choices=['assembly', 'reads', 'both'], default='assembly',
                        help="Validation source: 'assembly' (default, current behavior), 'reads' "
                             "(skip the assembly entirely and validate against BAM long reads), or "
                             "'both' (assembly first, reads to rescue misses).")
    parser.add_argument('--read_error_threshold', type=float, default=0.1,
                        help='Max edlib error rate for a read to validate a reconstruction (read modes). '
                             'Keep <= the assembly error threshold (0.1) for consistent hit/miss calls.')
    parser.add_argument('--min_read_support', type=int, default=1,
                        help='Min number of spanning reads that must clear --read_error_threshold (read modes).')
    parser.add_argument('--max_reads_per_site', type=int, default=1000,
                        help='Cap on candidate reads gathered per SV locus (read modes); '
                             'bounds work on deep read pileups.')
    parser.add_argument('--report', choices=['none', 'json'], default='none',
                        help="Write a per-SV evaluation sidecar. 'json' emits "
                             "<logbasename>.eval.jsonl next to the log (one JSON record per SV: "
                             "svid, outcome, tier, and per-segment status/reason/source/error/"
                             "junctions). Default 'none' (log and score table are unchanged).")
    args = parser.parse_args()

    os.makedirs('./logs', exist_ok=True)
    log_basename = log_name_from_config(args.config) or f'reconstruction_{timestamp}'
    # Keep eval modes in separate logs so a read-based run doesn't overwrite the
    # assembly baseline for the same experiment ('assembly' keeps the bare name).
    if args.eval_mode != 'assembly':
        log_basename += f'.{args.eval_mode}'
    log_filename = os.path.join('./logs', f'{log_basename}.log')
    report_path = os.path.join('./logs', f'{log_basename}.eval.jsonl') if args.report == 'json' else None

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

    logger.info(f'Logging to {log_filename}')
    logger.info(f'Config: {vars(args)}')

    if args.config:
        logger.info(f'Inferring params from {args.config}')
        args = update_args_from_config(args)

    global global_ref
    global global_sample
    global global_aligners
    global global_eval_mode, global_bam_reader, global_read_error_threshold
    global global_min_read_support, global_max_reads_per_site

    global_eval_mode = args.eval_mode
    global_read_error_threshold = args.read_error_threshold
    global_min_read_support = args.min_read_support
    global_max_reads_per_site = args.max_reads_per_site

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

    logger.info('Loading reference bytearrays')
    global_ref = load_fasta_to_bytes(args.reference, chroms)
    logger.info(f'Loaded {len(global_ref)} reference chromosomes.')

    global_sample = []
    global_aligners = defaultdict(list)

    if global_eval_mode in ('assembly', 'both'):
        logger.info('Loading sample assembly bytearrays')
        global_sample = [load_fasta_to_bytes(samp, chroms) for samp in args.sample]
        logger.info(f'Initialized {len(global_sample)} sample assembly file(s).')

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
    else:
        logger.info('Read-only eval mode: skipping sample assembly load and aligner build.')

    if global_eval_mode in ('reads', 'both'):
        if not args.bam:
            logger.error('--eval_mode reads/both requires a BAM (--bam) for read-based evaluation.')
            sys.exit(1)
        logger.info(f'Initializing thread-safe BAM reader: {args.bam}')
        global_bam_reader = BamReader(args.bam)

    precision, correct_calls, total_calls, inconclusive_calls, assembly_hits, read_hits = scorer.score_all(
        location_tolerance=args.location_tolerance, report_path=report_path)

    df = pd.DataFrame({
        'correct_calls': correct_calls,
        'total_calls': total_calls,
        'inconclusive': inconclusive_calls,
        'precision': precision,
        'assembly_hits': assembly_hits,
        'read_hits': read_hits,
    })
    # count columns: missing SV-type keys (e.g. a type with 0 hits) -> 0, not NaN
    for col in ('correct_calls', 'total_calls', 'inconclusive', 'assembly_hits', 'read_hits'):
        df[col] = df[col].fillna(0).astype(int)
    df.sort_values('total_calls', ascending=False, inplace=True)

    # to_string() prints all columns (default repr truncates the middle ones)
    logger.info('Score table:\n' + df.to_string())

    export_igv_session(args.calls, args.bam, args.classified, timestamp, args.igv_prefix)


if __name__ == '__main__':
    main()