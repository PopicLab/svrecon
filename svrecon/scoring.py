"""Scoring engine: the AlignScorer orchestrator and per-SV scoring."""
from dataclasses import dataclass, field
import json
import logging
import traceback
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Union, Tuple, Optional

import mappy
import pysam
from intervaltree import IntervalTree
from pysam import VariantRecord
from tqdm import tqdm

from svrecon.align import (check_match, edlib_to_cigartuples,
                           run_edlib_fallback, validate_segments_from_cigar, SeqSegmentValidationResult)
from svrecon.constants import Outcome, SubseqReason, SubseqStatus
from svrecon.plot import plot_dot_plot
from svrecon.reads import run_read_edlib, ReadEdlibResult
from svrecon.reconstruct import simulate_subsequences, QueryReconSubsequence
from svrecon.util import clamp, reverse_complement, merge_intervals

logger = logging.getLogger(__name__)


MIN_EDLIB_QUERY = 5000
# Search-window radius for the assembly edlib fallback (a legacy path for short
# queries mappy missed). This is a COST bound on an unbounded O(n*m) edlib search,
# kept separate from --location_tolerance: the latter only filters mappy hits by
# genomic position and defaults to infinity (meaningless for unscaffolded, per-
# contig assemblies where a hit's r_st is contig-local, not a genomic coordinate).
EDLIB_FALLBACK_MAX_TOLERANCE = 10_000_000
JUNCTION_VALIDATION_WINDOW_MAX = 5000

@dataclass
class ScorerRecord:
    lowest_error: float = 1.0
    passed: bool = False
    validating_seq: Optional[str] = None
    junction_results: List[SeqSegmentValidationResult] = field(default_factory=list)
    source: Optional[str] = None
    best_strand_match: Optional[int] = None

    # Failure diagnostics -- only meaningful once a tier call did NOT pass; used to
    # classify WHY a subsequence failed (junction rejection vs. no reads at all vs.
    # reads present but not spanning vs. spanning reads all aborted).
    seq_no_reads: bool = False
    seq_inconclusive: bool = False
    seq_reads_aborted: bool = False
    seq_had_junction_rejection: bool = False
    seq_best_fail_err: Optional[float] = None

    def update(self, **kwargs) -> None:
        if kwargs.get('passed'):
            if kwargs.get('lowest_error', self.lowest_error) >= self.lowest_error:
                return
            self.passed = True
            self.lowest_error = kwargs['lowest_error']
            self.validating_seq = kwargs.get('validating_seq')
            self.junction_results = kwargs.get('junction_results') or []
            self.source = kwargs.get('source')
            self.best_strand_match = kwargs.get('best_strand_match')
        else:
            if kwargs.get('no_reads'):
                self.seq_no_reads = True
            if kwargs.get('inconclusive'):
                self.seq_inconclusive = True
            if kwargs.get('reads_aborted'):
                self.seq_reads_aborted = True
            if kwargs.get('had_junction_rejection'):
                self.seq_had_junction_rejection = True
            best_fail_err = kwargs.get('best_fail_err')
            if best_fail_err is not None:
                self.seq_best_fail_err = best_fail_err if self.seq_best_fail_err is None else min(self.seq_best_fail_err, best_fail_err)


@dataclass
class AmbiguityResult:
    reference_match: bool = False
    reference_err: Optional[float] = None
    reference_strand: Optional[int] = None

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

        self.buffer = int(buffer)
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
        self.junction_window_factor = 0.05 # TODO make these constants config params
        self.junction_window_min = 150
        self.junction_window_max = JUNCTION_VALIDATION_WINDOW_MAX

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

    def score_all(self, location_tolerance=float('inf'), error_threshold=0.1, n_threads=40, report_path=None,
                 plot_first_n=0, plot_out_dir=None):
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

        # Optional dot-plot PNGs (recon subsequence vs. its validating real-data sequence),
        # up to plot_first_n per SV type. Off by default (plot_first_n=0).
        if plot_first_n:
            Path(plot_out_dir).mkdir(parents=True, exist_ok=True)

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
                    outcome = result.get('outcome', Outcome.HIT if result['is_correct'] else Outcome.MISS)
                    match_scores = result['match_scores']
                    line = f'{svid}\t{sv_type}\t{recon_sequences[0].chrom}:{coords[0]}-{coords[1]}\t{outcome.value}\t{match_scores}'
                    if outcome == Outcome.HIT:
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

                    if plot_first_n and total_calls[sv_type] <= plot_first_n:
                        # TODO: extend to cover every subsequence of a compound SV
                        validating_sequences = result.get('validating_sequences', [])
                        segments = result.get('segments', [])
                        first_recon_seq = recon_sequences[0].sequence if recon_sequences else None
                        first_validating_seq = validating_sequences[0] if validating_sequences else None
                        source = segments[0]['source'] if segments else None
                        if first_recon_seq is not None and first_validating_seq is not None:
                            plot_path = Path(plot_out_dir) / f'{sv_type}_{svid}_{source}.png'

                            logger.info(f'Producing dotplot for SV {svid} ({sv_type}, source={source})')
                            plot_dot_plot(first_recon_seq, first_validating_seq,
                                         f'{sv_type} {svid}', str(plot_path),
                                         s1_name='reconstructed', s2_name=f'validating ({source})')

                        else:
                            logger.warning(f'Skipping dot plot for {svid} ({sv_type}): no validating sequence available.')

                    if outcome == Outcome.HIT:
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
                    elif outcome == Outcome.INCONCLUSIVE:
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
                 match_error_threshold=0.1): 
        """Score one SV: reconstruct its alt allele(s) and validate each against reads/assembly
        per `self.eval_mode`, then roll up one outcome for the call ('hit' iff every reconstructed
        subsequence passed, 'miss' if any was contradicted, else 'inconclusive').

        Returns a dict with svid, sv_type, outcome, is_correct, tier, validation_source,
        rescued_by_edlib, coords, match_scores, validating_sequences, recon_sequences,
        subseq_diagnostics, and segments (the last for the optional `--report json` sidecar).
        """
        try:
            svid = records[0].info['SVID']
            sv_type = records[0].info['SVTYPE']

            recon_sequences: List[QueryReconSubsequence] = simulate_subsequences(records, buffer, self.ref)
            match_scores: List[float] = []
            validating_sequences: List[str] = []
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
            eval_reads = self.eval_mode in ('reads', 'both')
            eval_assembly = self.eval_mode in ('assembly', 'both')

            
            # if eval_reads: # NOTE: why cant we take the min / max of sv_ref_span_pts?
            #     starts = [r.start for r in records]
            #     stops = [r.stop for r in records]
            #     for r in records:
            #         if 'TARGET' in r.info:
            #             starts.append(int(r.info['TARGET']))
            #             stops.append(int(r.info['TARGET']))
            #     sv_ref_start_bp = max(0, min(starts))
            #     sv_ref_end_bp = max(stops)
            candidate_reads: Dict = self.bam_reader.candidate_reads(records) if eval_reads else None # TODO: 

            for query in recon_sequences:
                # Window size 
                junction_validation_radius = int(clamp(round(self.junction_window_factor * len(query.sequence)), self.junction_window_min, self.junction_window_max))

                score_record = ScorerRecord()

                best_score = 1.0
                mappy_matched = False

                seq_had_junction_rejection = False
                seq_best_rejected_err = 1.0

                # seq_has_aligner = chrom in self.aligners
                # seq_saw_candidate = False
                # seq_best_fail_err = None
                # seq_no_reads = False
                # seq_inconclusive = False  # reads overlap but none span the full resulting allele
                # seq_reads_aborted = False  # spanning reads existed but all exceeded the 2x edlib bound
                # seq_source = None  # which tier validated this subsequence: assembly | edlib | reads
                # seq_passed = False  # a tier validated it within ITS OWN threshold (+ junctions)
                # seq_junctions = []  # per-junction results from the tier whose check we report
                # seq_strand = None  # strand of the assembly match (meaningful only there; None for reads/edlib)
                # seq_validating_sequence = None  # the real-data sequence (read/assembly span/edlib span) that validated this subsequence

                # Read-based evaluation
           

                if eval_reads: 
                    read_based_score = self._score_read_based_eval(query, junction_validation_radius, sv_ref_start_bp, sv_ref_end_bp, buffer)
                    score_record.update(**read_based_score)

                if not score_record.passed and eval_assembly:
                    assembly_based_score = self._score_assembly_based_eval(query, junction_validation_radius, location_tolerance, match_error_threshold)
                    score_record.update(**assembly_based_score)

                    if not score_record.passed and len(query.sequence) < MIN_EDLIB_QUERY:
                        edlib_based_score = self._score_edlib_based_eval(query, junction_validation_radius, match_error_threshold, buffer)
                        score_record.update(**edlib_based_score)

                ambiguity_check_results = None
                if score_record.passed and self.check_reference:
                    ambiguity_check_results = self._check_ambiguity(query, junction_validation_radius, match_error_threshold, buffer)
                    

                # if self.eval_mode in ('reads', 'both'):
                #     read_seqs = self.bam_reader.candidate_read_seqs(
                #         chrom, sv_ref_start_bp, sv_ref_end_bp, int(buffer), self.max_reads_per_site)
                #     # A read can only confirm the variant if it is long enough to hold
                #     # the resulting allele
                #     result_len = query.result_len or len(sequence) # NOTE
                #     spanning = [r for r in read_seqs if len(r) >= result_len]
                #     if not read_seqs:
                #         seq_no_reads = True
                #     elif not spanning:
                #         seq_inconclusive = True
                #     else:
                #         read_res = run_read_edlib(sequence, spanning,
                #                                   self.read_error_threshold, self.min_read_support)
                #         if read_res.cigar is not None: # If there is a valid read
                #             seq_saw_candidate = True
                #             if read_res.error <= self.read_error_threshold:
                #                 cigartuples = edlib_to_cigartuples(read_res.cigar)
                #                 junctions_validation_results: List[SeqSegmentValidationResult] = \
                #                     validate_segments_from_cigar(cigartuples, query.segments,
                #                                                  radius=junction_window_size,
                #                                                  error_threshold=self.read_error_threshold)
                #                 if not seq_passed:
                #                     seq_junctions = junctions_validation_results
                #                 if all(j.passed for j in junctions_validation_results):
                #                     best_score = min(best_score, read_res.error)
                #                     seq_source = 'reads'
                #                     seq_passed = True
                #                     seq_validating_sequence = read_res.matched_read_sequence
                #                 else:
                #                     seq_had_junction_rejection = True
                #                     seq_best_rejected_err = min(seq_best_rejected_err, read_res.error)
                #                     seq_best_fail_err = read_res.error if seq_best_fail_err is None else min(seq_best_fail_err, read_res.error)
                #             else:
                #                 seq_best_fail_err = read_res.error if seq_best_fail_err is None else min(seq_best_fail_err, read_res.error)
                #         else:
                #             # spanning reads existed but all exceeded the 2x edlib bound
                #             seq_reads_aborted = True

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

                        # forward strand only user
                        for a in alignments:
                            seq_saw_candidate = True
                            err = check_match(a, sequence)
                            if err <= match_error_threshold:
                                junctions_validation_results: List[SeqSegmentValidationResult] = \
                                    validate_segments_from_cigar(a.cigar, query.segments,
                                    radius=junction_window_size,
                                    q_st=a.q_st, q_en=a.q_en,
                                    query_len=len(sequence), strand=a.strand,
                                    error_threshold=match_error_threshold)
                                
                                if not seq_passed: # NOTE: why this checked?
                                    seq_junctions = junctions_validation_results
                                if all(j.passed for j in junctions_validation_results):
                                    if err < best_score:
                                        best_score = err
                                        seq_strand = a.strand   # record the best passing match's strand
                                        matched = aligner.seq(a.ctg, a.r_st, a.r_en)
                                        # seq_validating_sequence = matched # Assume forward strand to forward strand match
                                        if a.strand == -1: logger.warning(f"{svid} {sv_type} mapped to reverse strand") # We shouldn't be mappign via reverse?
                                        seq_validating_sequence = reverse_complement(matched) if a.strand == -1 else matched
                                    mappy_matched = True
                                    seq_source = 'assembly'
                                    seq_passed = True
                                else:
                                    seq_had_junction_rejection = True
                                    seq_best_rejected_err = min(seq_best_rejected_err, err)
                                    seq_best_fail_err = err if seq_best_fail_err is None else min(seq_best_fail_err, err)
                            else:
                                seq_best_fail_err = err if seq_best_fail_err is None else min(seq_best_fail_err, err)

                # if self.eval_mode in ('assembly', 'both') and not seq_passed and len(sequence) < MIN_EDLIB_QUERY and not mappy_matched:
                #     # Assembly edlib fallback: for short subsequences that mappy
                #     # failed to align (mappy can miss very short queries), retry with
                #     # edlib against expanding windows of the assembly, applying the
                #     # same error-threshold + junction-validation acceptance criteria.
                #     mappy_passed_all = False
                #     for samp_bytes in self.sample:
                #         edlib_res = run_edlib_fallback(
                #             sequence,
                #             chrom,
                #             query.location,
                #             int(buffer),
                #             EDLIB_FALLBACK_MAX_TOLERANCE,
                #             match_error_threshold,
                #             samp_bytes
                #         )
                #         if edlib_res:
                #             seq_saw_candidate = True
                #             if edlib_res.error <= match_error_threshold:
                #                 cigartuples = edlib_to_cigartuples(edlib_res.cigar)
                #                 junctions_validation_results: List[SeqSegmentValidationResult] = \
                #                     validate_segments_from_cigar(cigartuples, query.segments,
                #                                                  radius=junction_window_size,
                #                                                  error_threshold=match_error_threshold)
                #                 if not seq_passed:
                #                     seq_junctions = junctions_validation_results
                #                 if all(j.passed for j in junctions_validation_results):
                #                     if edlib_res.error < best_score:
                #                         best_score = edlib_res.error
                #                         seq_validating_sequence = edlib_res.matched_target_sequence
                #                     seq_source = 'edlib'
                #                     seq_passed = True
                #                 else:
                #                     seq_had_junction_rejection = True
                #                     seq_best_rejected_err = min(seq_best_rejected_err, edlib_res.error)
                #                     seq_best_fail_err = edlib_res.error if seq_best_fail_err is None else min(seq_best_fail_err, edlib_res.error)
                #             else:
                #                 seq_best_fail_err = edlib_res.error if seq_best_fail_err is None else min(seq_best_fail_err, edlib_res.error)


                match_scores.append(score_record.lowest_error)
                validating_sequences.append(score_record.validating_seq)

                # Optional ambiguity check: does the reconstructed allele ALSO validate against the
                # REFERENCE (which lacks the SV) under the SAME requirements as the assembly check --
                # bulk error AND junction validation (flanking context)? If so, the SV's full
                # structure exists in the reference too, so the assembly/read hit is not specific to
                # the SV (common in repetitive / segmental-dup / inverted-repeat regions) -> the call
                # is ambiguous, not a hit. mappy is both-strand + genome-wide; edlib is the
                # short-allele fallback only (and forward-only). Records the reference match strand.
                # seq_reference_match = False
                # seq_reference_err = None
                # seq_reference_strand = None
                # if seq_passed and self.check_reference:
                #     if self.reference_aligners:
                #         for aligner in self.reference_aligners.get(chrom, []):
                #             for a in aligner.map(sequence):
                #                 err = check_match(a, sequence)
                #                 if err <= match_error_threshold and (seq_reference_err is None or err < seq_reference_err):
                #                     junctions_validation_results: List[SeqSegmentValidationResult] = \
                #                         validate_segments_from_cigar(a.cigar, query.segments,
                #                                                      radius=junction_window_size,
                #                                                      q_st=a.q_st, q_en=a.q_en,
                #                                                      query_len=len(sequence), strand=a.strand,
                #                                                      error_threshold=match_error_threshold)
                #                     if all(j.passed for j in junctions_validation_results):
                #                         seq_reference_match = True
                #                         seq_reference_err = err
                #                         seq_reference_strand = a.strand
                #     if not seq_reference_match and len(sequence) < MIN_EDLIB_QUERY:
                #         ref_res = run_edlib_fallback(sequence, chrom, query.location, int(buffer),
                #                                      self.reference_search_tolerance, match_error_threshold, self.ref)
                #         if ref_res is not None and ref_res.error <= match_error_threshold:
                #             junctions_validation_results: List[SeqSegmentValidationResult] = \
                #                 validate_segments_from_cigar(edlib_to_cigartuples(ref_res.cigar),
                #                                              query.segments,
                #                                              radius=junction_window_size,
                #                                              error_threshold=match_error_threshold)
                #             if all(j.passed for j in junctions_validation_results):
                #                 seq_reference_match = True
                #                 seq_reference_err = ref_res.error
                #                 seq_reference_strand = 1  # edlib aligns forward only

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
                status = reason = None

                # TODO: refactor so the result log just shows all the flags to avoid having code to interpret it

                if score_record.passed and ambiguity_check_results is not None and ambiguity_check_results.reference_match:
                    status, reason = SubseqStatus.INCONCLUSIVE, SubseqReason.REFERENCE_MATCH
                elif score_record.passed:
                    status, reason = SubseqStatus.PASS, SubseqReason.PASS
                elif score_record.seq_had_junction_rejection or seq_saw_candidate or score_record.seq_reads_aborted:
                    status = SubseqStatus.FAIL
                    reason = SubseqReason.JUNCTION_FAILED if score_record.seq_had_junction_rejection else SubseqReason.OVER_ERROR_THRESHOLD
                elif score_record.seq_no_reads or score_record.seq_inconclusive:
                    status = SubseqStatus.INCONCLUSIVE
                    reason = SubseqReason.NO_READS if score_record.seq_no_reads else SubseqReason.INCONCLUSIVE
                else:
                    status = SubseqStatus.FAIL
                    reason = SubseqReason.NO_ALIGNER if not seq_has_aligner else SubseqReason.NO_ALIGNMENT_IN_WINDOW

                if reason == SubseqReason.REFERENCE_MATCH:
                    err_str = f'{ambiguity_check_results.reference_err:.4f}'  # error of the coincidental reference match
                elif score_record.seq_best_fail_err is not None:
                    err_str = f'{score_record.seq_best_fail_err:.4f}'
                elif reason == SubseqReason.OVER_ERROR_THRESHOLD:
                    err_str = 'aborted'  # spanning reads exceeded the 2x-threshold edlib bound
                else:
                    err_str = 'NA'
                subseq_diagnostics.append(f'{reason.value}:{err_str}')
                subseq_sources.append(seq_source)
                subseq_status.append(status)

                # Structured per-subsequence record for the optional JSON report. `error` is the
                # winning error when the subsequence passed, else the best failing error, else None
                # (nothing testable produced a number).
                if status == SubseqStatus.PASS:
                    seg_error = float(f'{best_score:.4g}')
                elif reason == SubseqReason.REFERENCE_MATCH:
                    seg_error = float(f'{seq_reference_err:.4g}')
                elif seq_best_fail_err is not None:
                    seg_error = float(f'{seq_best_fail_err:.4g}')
                else:
                    seg_error = None
                segments.append({
                    'chrom': chrom,
                    'ref_start': query.location,
                    'status': status.value,
                    'reason': reason.value,
                    'source': seq_source,
                    'error': seg_error,
                    # Strand of the decisive alignment where meaningful: the reference match for a
                    # reference_match, else the assembly match. null for read/edlib confirmations
                    # (reads are randomly oriented; strand carries no signal there).
                    'strand': seq_reference_strand if reason == SubseqReason.REFERENCE_MATCH else seq_strand,
                    'junctions': [j.jsonify() for j in seq_junctions],
                    'validating_sequence': seq_validating_sequence,
                })

                # TODO: what is this comparison doing? why is it > match_error_threshold?
                if score_record.lowest_error > match_error_threshold and score_record.seq_had_junction_rejection:
                    sv_junction_rejected = True
                    best_rejected_err = score_record.seq_best_fail_err

            coords = records[0].pos, records[0].stop

            # SV-level roll-up:
            #   hit          -> every subsequence passed
            #   miss         -> at least one subsequence was contradicted (spanning
            #                   reads / assembly aligned but didn't match)
            #   inconclusive -> no contradiction, but at least one subsequence could
            #                   not be tested (reads mode: no read spans the full allele)
            if len(subseq_status) > 0 and all(s == SubseqStatus.PASS for s in subseq_status):
                outcome = Outcome.HIT
            elif SubseqStatus.FAIL in subseq_status:
                outcome = Outcome.MISS
            elif SubseqStatus.INCONCLUSIVE in subseq_status:
                outcome = Outcome.INCONCLUSIVE
            else:
                outcome = Outcome.MISS
            is_correct = (outcome == Outcome.HIT)

            if outcome == Outcome.MISS and sv_junction_rejected:
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
                'validating_sequences': validating_sequences,
                'recon_sequences': recon_sequences,
                'subseq_diagnostics': subseq_diagnostics,
                'segments': segments,
            }
        except Exception as e:
            error_msg = f'Worker failed for SVID: {records[0].info.get("SVID", "Unknown")}. Error: {e}\n{traceback.format_exc()}'
            logger.error(error_msg)
            raise e
        
    def _score_read_based_eval(self, query: QueryReconSubsequence, junction_validation_radius, read_bp_start, read_bp_end, buffer) -> Tuple[Dict, bool]:
        '''
        Read based validation,
        '''
        lowest_error = 1.0
        passed = False
        validating_seq = None
        junctions_validation_results = None

        no_reads = False
        inconclusive = False
        reads_aborted = False
        had_junction_rejection = False
        best_fail_err = None

        reads_overlapping = self.bam_reader.candidate_read_seqs(query.chrom, read_bp_start, read_bp_end, int(buffer), self.max_reads_per_site)
        reads_spanning = [r for r in reads_overlapping if len(r) >= (query.result_len or len(query.sequence))] # TODO: I think this is just query.sequence

        if not reads_overlapping:
            no_reads = True
        elif not reads_spanning:
            inconclusive = True
        else:
            read_res: ReadEdlibResult = run_read_edlib(query.sequence, reads_spanning, self.read_error_threshold, self.min_read_support)

            if read_res.is_read_found():
                if read_res.error <= self.read_error_threshold:
                    cigartuples = edlib_to_cigartuples(read_res.cigar)
                    junctions_validation_results: List[SeqSegmentValidationResult] = \
                                                                        validate_segments_from_cigar(cigartuples, query.segments,
                                                                                                radius=junction_validation_radius,
                                                                                                error_threshold=self.read_error_threshold)
                    if all(j.passed for j in junctions_validation_results):
                        validating_seq = read_res.matched_read_sequence
                        passed = True
                        lowest_error = read_res.error
                    else:
                        had_junction_rejection = True
                        best_fail_err = read_res.error
                else:
                    best_fail_err = read_res.error
            else:
                # spanning reads existed but all exceeded the 2x edlib bound
                reads_aborted = True
        return {
            'source': 'reads',
            'lowest_error': lowest_error,
            'passed': passed,
            'validating_seq': validating_seq,
            'junction_results': junctions_validation_results,
            'no_reads': no_reads,
            'inconclusive': inconclusive,
            'reads_aborted': reads_aborted,
            'had_junction_rejection': had_junction_rejection,
            'best_fail_err': best_fail_err,
        }

    def _score_assembly_based_eval(self, query: QueryReconSubsequence, junction_validation_radius, location_tolerance, match_error_threshold):
        lowest_error = 1.0
        passed = False
        validating_seq = None
        junctions_validation_results = None
        had_junction_rejection = False
        best_fail_err = None
        # Assembly match specific fields
        best_strand_match = None

        for aligner in self.aligners[query.chrom]:
            all_alignments: List[mappy.Alignment] = list(aligner.map(query.sequence))
            alignments = [a for a in all_alignments if abs(a.r_st - query.location) <= location_tolerance]

            for alignment in alignments:
                match_err = check_match(alignment, query.sequence)
                if match_err <= match_error_threshold:
                    junctions_validation_results: List[SeqSegmentValidationResult] = \
                                                        validate_segments_from_cigar(alignment.cigar, query.segments,
                                                        radius=junction_validation_radius,
                                                        q_st=alignment.q_st, q_en=alignment.q_en,
                                                        query_len=len(query.sequence), strand=alignment.strand,
                                                        error_threshold=match_error_threshold)
                    if all(j.passed for j in junctions_validation_results):
                        if match_err < lowest_error:
                            passed = True
                            lowest_error = match_err
                            best_strand_match = alignment.strand   # record the best passing match's strand
                            validating_seq = aligner.seq(alignment.ctg, alignment.r_st, alignment.r_en)
                            if alignment.strand == -1: logger.warning(f"mapped to reverse strand") # NOTE: We shouldn't be mappign via reverse?
                            # validating_seq = reverse_complement(best_seq_match) if a.strand == -1 else best_seq_match
                    else:
                        had_junction_rejection = True
                        best_fail_err = match_err if best_fail_err is None else min(best_fail_err, match_err)
                else:
                    best_fail_err = match_err if best_fail_err is None else min(best_fail_err, match_err)
        return {
            'source': 'assembly',
            'lowest_error': lowest_error,
            'passed': passed,
            'validating_seq': validating_seq,
            'junction_results': junctions_validation_results,
            # assembly match specific fields
            'best_strand_match': best_strand_match,
            'had_junction_rejection': had_junction_rejection,
            'best_fail_err': best_fail_err,
        }

    def _score_edlib_based_eval(self, query: QueryReconSubsequence, junction_validation_radius, match_error_threshold, buffer):
        lowest_error = 1.0
        passed = False
        validating_seq = None
        junctions_validation_results = None

        had_junction_rejection = False
        best_fail_err = None

        for samp_bytes in self.sample:
            edlib_res = run_edlib_fallback(
                query.sequence,
                query.chrom,
                query.location,
                buffer,
                EDLIB_FALLBACK_MAX_TOLERANCE,
                match_error_threshold,
                samp_bytes
            )
            if edlib_res and edlib_res.error <= match_error_threshold:
                cigartuples = edlib_to_cigartuples(edlib_res.cigar)
                junctions_validation_results: List[SeqSegmentValidationResult] = \
                    validate_segments_from_cigar(cigartuples, query.segments,
                                                    radius=junction_validation_radius,
                                                    error_threshold=match_error_threshold)
                if all(j.passed for j in junctions_validation_results):
                    if edlib_res.error < lowest_error:
                        lowest_error = edlib_res.error
                        validating_seq = edlib_res.matched_target_sequence
                        passed = True
                else:
                    had_junction_rejection = True
                    best_fail_err = edlib_res.error if best_fail_err is None else min(best_fail_err, edlib_res.error)
            elif edlib_res:
                best_fail_err = edlib_res.error if best_fail_err is None else min(best_fail_err, edlib_res.error)

        return {
            'source': 'edlib',
            'lowest_error': lowest_error,
            'passed': passed,
            'validating_seq': validating_seq,
            'junction_results': junctions_validation_results,
            'had_junction_rejection': had_junction_rejection,
            'best_fail_err': best_fail_err,
        }

    def _check_ambiguity(self, query, junction_window_size, match_error_threshold, buffer):
        seq_reference_match = False
        seq_reference_err = None
        seq_reference_strand = None
        if self.reference_aligners:
            for aligner in self.reference_aligners.get(query.chrom, []):
                for a in aligner.map(query.sequence):
                    err = check_match(a, query.sequence)
                    if err <= match_error_threshold and (seq_reference_err is None or err < seq_reference_err):
                        junctions_validation_results: List[SeqSegmentValidationResult] = \
                            validate_segments_from_cigar(a.cigar, query.segments,
                                                            radius=junction_window_size,
                                                            q_st=a.q_st, q_en=a.q_en,
                                                            query_len=len(query.sequence), strand=a.strand,
                                                            error_threshold=match_error_threshold)
                        if all(j.passed for j in junctions_validation_results):
                            seq_reference_match = True
                            seq_reference_err = err
                            seq_reference_strand = a.strand
        if not seq_reference_match and len(query.sequence) < MIN_EDLIB_QUERY:
            ref_res = run_edlib_fallback(query.sequence, query.chrom, query.location, int(buffer),
                                            self.reference_search_tolerance, match_error_threshold, self.ref)
            if ref_res is not None and ref_res.error <= match_error_threshold:
                junctions_validation_results: List[SeqSegmentValidationResult] = \
                    validate_segments_from_cigar(edlib_to_cigartuples(ref_res.cigar),
                                                    query.segments,
                                                    radius=junction_window_size,
                                                    error_threshold=match_error_threshold)
                if all(j.passed for j in junctions_validation_results):
                    seq_reference_match = True
                    seq_reference_err = ref_res.error
                    seq_reference_strand = 1  # edlib aligns forward only
        return AmbiguityResult(
            reference_match=seq_reference_match,
            reference_err=seq_reference_err,
            reference_strand=seq_reference_strand,
        )

    # TODO:
    # - track failed junctions