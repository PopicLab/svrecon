"""Scoring engine: the AlignScorer orchestrator and per-SV scoring."""
import json
import logging
import traceback
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Dict, List, Union

import mappy
import pysam
from intervaltree import IntervalTree
from pysam import VariantRecord
from tqdm import tqdm

from svrecon.align import (check_match, edlib_to_cigartuples,
                           run_edlib_fallback, validate_junctions_from_cigar, SeqJunctionsValidationResult)
from svrecon.reads import run_read_edlib
from svrecon.reconstruct import simulate_subsequences, QueryReconSubsequence

logger = logging.getLogger(__name__)


MIN_EDLIB_QUERY = 5000
# Search-window radius for the assembly edlib fallback (a legacy path for short
# queries mappy missed). This is a COST bound on an unbounded O(n*m) edlib search,
# kept separate from --location_tolerance: the latter only filters mappy hits by
# genomic position and defaults to infinity (meaningless for unscaffolded, per-
# contig assemblies where a hit's r_st is contig-local, not a genomic coordinate).
EDLIB_FALLBACK_MAX_TOLERANCE = 10_000_000
JUNCTION_VALIDATION_WINDOW = 300


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
        # Run-time state, set by the CLI after ref/aligners are built. Shared by
        # reference across worker threads (ThreadPoolExecutor), never copied.
        self.ref = None
        self.sample = []
        self.aligners = None
        self.eval_mode = 'assembly'
        self.bam_reader = None
        self.read_error_threshold = 0.1
        self.min_read_support = 1
        self.max_reads_per_site = 1000
        # Optional ambiguity check: also search the REFERENCE for the reconstructed allele;
        # a hit there means the match isn't specific to the SV -> inconclusive. Off by default.
        self.check_reference = False
        self.reference_aligners = None                 # mappy aligners over the reference
        self.reference_search_tolerance = EDLIB_FALLBACK_MAX_TOLERANCE
        # Junction-validation window, scaled to SV size and clamped: window = clamp(factor*span,
        # min, max). Applies to ALL junction checks (reads, assembly, edlib, reference) so the
        # junction test measures the SV's actual change rather than the surrounding context --
        # otherwise a small SV's signal is diluted below threshold and passes on flanks alone.
        self.junction_window_factor = 1.5
        self.junction_window_min = 150
        self.junction_window_max = JUNCTION_VALIDATION_WINDOW

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
                    self.score_sv,
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
                    recon_sequences = result['recon_sequences']
                    coords = result['coords']
                    outcome = result.get('outcome', 'hit' if result['is_correct'] else 'miss')
                    match_scores = result['match_scores']
                    line = f'{svid}\t{sv_type}\t{recon_sequences[0]["chrom"]}:{coords[0]}-{coords[1]}\t{outcome}\t{match_scores}'
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

    def score_sv(self, records: List[VariantRecord], buffer: Union[int, float], location_tolerance: Union[int, float],
                 error_threshold=0.1): 
        """Score one SV: reconstruct its alt allele(s) and validate each against reads/assembly
        per `self.eval_mode`, then roll up one outcome for the call ('hit' iff every reconstructed
        subsequence passed, 'miss' if any was contradicted, else 'inconclusive').

        Returns a dict with svid, sv_type, outcome, is_correct, tier, validation_source,
        rescued_by_edlib, coords, match_scores, recon_sequences, subseq_diagnostics, and segments
        (the last for the optional `--report json` sidecar).
        """
        try:
            svid = records[0].info['SVID']
            sv_type = records[0].info['SVTYPE']

            # Junction-validation window scaled to this SV's change magnitude (its reference
            # footprint incl. any TARGET) and clamped. Used by every junction check below so a
            # small SV's junction signal isn't diluted by the surrounding context window.
            span_pts = [p for r in records for p in (r.start, r.stop)] # NOTE: are span pts effectively duplicated by insilicoSV operations?
            span_pts += [int(r.info['TARGET']) for r in records if 'TARGET' in r.info]
            sv_span = (max(span_pts) - min(span_pts)) if span_pts else 0
            junction_window = int(min(self.junction_window_max,
                                      max(self.junction_window_min,
                                          round(self.junction_window_factor * sv_span))))

            recon_sequences: List[QueryReconSubsequence] = simulate_subsequences(records, buffer, self.ref)
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
            if self.eval_mode in ('reads', 'both'):
                starts = [r.start for r in records]
                stops = [r.stop for r in records]
                for r in records:
                    if 'TARGET' in r.info:
                        starts.append(int(r.info['TARGET']))
                        stops.append(int(r.info['TARGET']))
                read_bp_start = max(0, min(starts))
                read_bp_end = max(stops)

            for query in recon_sequences:
                sequence = query.sequence
                chrom = query.chrom

                best_score = 1.0
                mappy_matched = False

                seq_had_junction_rejection = False
                seq_best_rejected_err = 1.0

                seq_has_aligner = chrom in self.aligners
                seq_saw_candidate = False
                seq_best_fail_err = None
                seq_no_reads = False
                seq_inconclusive = False  # reads overlap but none span the full resulting allele
                seq_reads_aborted = False  # spanning reads existed but all exceeded the 2x edlib bound
                seq_source = None  # which tier validated this subsequence: assembly | edlib | reads
                seq_passed = False  # a tier validated it within ITS OWN threshold (+ junctions)
                seq_junctions = []  # per-junction results from the tier whose check we report
                seq_strand = None  # strand of the assembly match (meaningful only there; None for reads/edlib)

                # Read-based evaluation
                if self.eval_mode in ('reads', 'both'):
                    read_seqs = self.bam_reader.candidate_read_seqs(
                        chrom, read_bp_start, read_bp_end, int(buffer), self.max_reads_per_site)
                    # A read can only confirm the variant if it is long enough to hold
                    # the resulting allele
                    result_len = query.result_len or len(sequence)
                    spanning = [r for r in read_seqs if len(r) >= result_len]
                    if not read_seqs:
                        seq_no_reads = True
                    elif not spanning:
                        seq_inconclusive = True
                    else:
                        read_res, _ = run_read_edlib(sequence, spanning,
                                                     self.read_error_threshold, self.min_read_support)
                        if read_res:
                            seq_saw_candidate = True
                            if read_res['error'] <= self.read_error_threshold:
                                cigartuples = edlib_to_cigartuples(read_res['cigar'])
                                junctions_validation_results = validate_junctions_from_cigar(cigartuples, query.junctions,
                                                                 window=junction_window,
                                                                 error_threshold=self.read_error_threshold)
                                if not seq_passed:
                                    seq_junctions = junctions_validation_results
                                if all(j.passed for j in junctions_validation_results):
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

                if self.eval_mode in ('assembly', 'both') and not seq_passed and chrom in self.aligners:
                    # Assembly check (Tier 2). In 'both' mode this only runs when the
                    # reads above did not already validate the subsequence. Align the
                    # reconstructed subsequence to the per-chromosome assembly
                    # aligner(s), keep alignments within location_tolerance of the
                    # expected locus, and accept a match only when the error rate
                    # clears error_threshold AND the SV junctions validate.
                    for aligner in self.aligners[chrom]:
                        all_alignments: List[mappy.Alignment] = list(aligner.map(sequence))
                        alignments = [a for a in all_alignments if abs(a.r_st - query.location) <= location_tolerance]

                        for a in alignments:
                            seq_saw_candidate = True
                            err = check_match(a, sequence)
                            if err <= error_threshold:
                                junctions_validation_results: List[SeqJunctionsValidationResult] = \
                                    validate_junctions_from_cigar(a.cigar, query.junctions,
                                    window=junction_window,
                                    q_st=a.q_st, q_en=a.q_en,
                                    query_len=len(sequence), strand=a.strand,
                                    error_threshold=error_threshold)
                                
                                if not seq_passed: # NOTE: why this checked?
                                    seq_junctions = junctions_validation_results
                                if all(j.passed for j in junctions_validation_results):
                                    if err < best_score:
                                        best_score = err
                                        seq_strand = a.strand   # record the best passing match's strand
                                    mappy_matched = True
                                    seq_source = 'assembly'
                                    seq_passed = True
                                else:
                                    seq_had_junction_rejection = True
                                    seq_best_rejected_err = min(seq_best_rejected_err, err)
                                    seq_best_fail_err = err if seq_best_fail_err is None else min(seq_best_fail_err, err)
                            else:
                                seq_best_fail_err = err if seq_best_fail_err is None else min(seq_best_fail_err, err)

                if self.eval_mode in ('assembly', 'both') and not seq_passed and len(sequence) < MIN_EDLIB_QUERY and not mappy_matched:
                    # Assembly edlib fallback: for short subsequences that mappy
                    # failed to align (mappy can miss very short queries), retry with
                    # edlib against expanding windows of the assembly, applying the
                    # same error-threshold + junction-validation acceptance criteria.
                    mappy_passed_all = False
                    for samp_bytes in self.sample:
                        edlib_res = run_edlib_fallback(
                            sequence,
                            chrom,
                            query.location,
                            int(buffer),
                            EDLIB_FALLBACK_MAX_TOLERANCE,
                            error_threshold,
                            samp_bytes
                        )
                        if edlib_res:
                            seq_saw_candidate = True
                            if edlib_res['error'] <= error_threshold:
                                cigartuples = edlib_to_cigartuples(edlib_res['cigar'])
                                junctions_validation_results = validate_junctions_from_cigar(cigartuples, query.junctions,
                                                                 window=junction_window,
                                                                 error_threshold=error_threshold)
                                if not seq_passed:
                                    seq_junctions = junctions_validation_results
                                if all(j.passed for j in junctions_validation_results):
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

                # Optional ambiguity check: does the reconstructed allele ALSO validate against the
                # REFERENCE (which lacks the SV) under the SAME requirements as the assembly check --
                # bulk error AND junction validation (flanking context)? If so, the SV's full
                # structure exists in the reference too, so the assembly/read hit is not specific to
                # the SV (common in repetitive / segmental-dup / inverted-repeat regions) -> the call
                # is ambiguous, not a hit. mappy is both-strand + genome-wide; edlib is the
                # short-allele fallback only (and forward-only). Records the reference match strand.
                seq_reference_match = False
                seq_reference_err = None
                seq_reference_strand = None
                if seq_passed and self.check_reference:
                    if self.reference_aligners:
                        for aligner in self.reference_aligners.get(chrom, []):
                            for a in aligner.map(sequence):
                                err = check_match(a, sequence)
                                if err <= error_threshold and (seq_reference_err is None or err < seq_reference_err):
                                    junctions_validation_results = validate_junctions_from_cigar(a.cigar, query.junctions,
                                                                     window=junction_window,
                                                                     q_st=a.q_st, q_en=a.q_en,
                                                                     query_len=len(sequence), strand=a.strand,
                                                                     error_threshold=error_threshold)
                                    if all(j.passed for j in junctions_validation_results):
                                        seq_reference_match = True
                                        seq_reference_err = err
                                        seq_reference_strand = a.strand
                    if not seq_reference_match and len(sequence) < MIN_EDLIB_QUERY:
                        ref_res = run_edlib_fallback(sequence, chrom, query.location, int(buffer),
                                                     self.reference_search_tolerance, error_threshold, self.ref)
                        if ref_res is not None and ref_res['error'] <= error_threshold:
                            junctions_validation_results = validate_junctions_from_cigar(edlib_to_cigartuples(ref_res['cigar']),
                                                             query.junctions,
                                                             window=junction_window,
                                                             error_threshold=error_threshold)
                            if all(j.passed for j in junctions_validation_results):
                                seq_reference_match = True
                                seq_reference_err = ref_res['error']
                                seq_reference_strand = 1  # edlib aligns forward only

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
                if seq_passed and seq_reference_match:
                    status, reason = 'inconclusive', 'reference_match'
                elif seq_passed:
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

                if reason == 'reference_match':
                    err_str = f'{seq_reference_err:.4f}'  # error of the coincidental reference match
                elif seq_best_fail_err is not None:
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
                elif reason == 'reference_match':
                    seg_error = float(f'{seq_reference_err:.4g}')
                elif seq_best_fail_err is not None:
                    seg_error = float(f'{seq_best_fail_err:.4g}')
                else:
                    seg_error = None
                segments.append({
                    'chrom': chrom,
                    'ref_start': query.location,
                    'status': status,
                    'reason': reason,
                    'source': seq_source,
                    'error': seg_error,
                    # Strand of the decisive alignment where meaningful: the reference match for a
                    # reference_match, else the assembly match. null for read/edlib confirmations
                    # (reads are randomly oriented; strand carries no signal there).
                    'strand': seq_reference_strand if reason == 'reference_match' else seq_strand,
                    'junctions': [j.jsonify() for j in seq_junctions],
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
                'recon_sequences': recon_sequences,
                'subseq_diagnostics': subseq_diagnostics,
                'segments': segments,
            }
        except Exception as e:
            error_msg = f'Worker failed for SVID: {records[0].info.get("SVID", "Unknown")}. Error: {e}\n{traceback.format_exc()}'
            logger.error(error_msg)
            raise e
