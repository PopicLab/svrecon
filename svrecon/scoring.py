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
from svrecon.constants import *
from svrecon.plot import plot_query_dot_plots
from svrecon.reads import run_read_edlib, ReadEdlibResult
from svrecon.reconstruct import simulate_subsequences, QueryReconSubsequence
from svrecon.util import clamp, reverse_complement, get_start_stop

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
class QueryInfo:
    query: QueryReconSubsequence

    # Pass Diagnostics
    lowest_pass_error: float = 1.0
    lowest_error: float = 1.0  # lowest error seen overall, pass or fail
    passed: bool = False
    validating_seq: Optional[str] = None
    ref_sequence: Optional[str] = None
    
    junction_results: List[SeqSegmentValidationResult] = field(default_factory=list)
    source: Optional[ValidationSource] = None
    # Pass Diagnostics, Optional
    best_strand_match: Optional[int] = None

    # Failure diagnostics -- only meaningful once a tier call did NOT pass. Only the
    # reads tier emits these; default True ("not applicable/assume satisfied") so an
    # assembly/edlib-only call never spuriously trips the classification below.
    overlapping_spanning_reads_found: bool = True
    candidate_read_found: bool = True

    # Populated by apply_ambiguity() -- only meaningful when passed, and only checked
    # when --check_reference is on.
    reference_check_match: bool = False
    reference_check_err: Optional[float] = None
    reference_check_strand: Optional[int] = None

    # Populated by get_summary() -- the final per-subsequence verdict.
    status: Optional['SubseqStatus'] = None
    reason: Optional['SubseqReason'] = None

    def update(self, **kwargs) -> None:
        # Field update based on pass / no pass on validation methods
        self.lowest_error = min(self.lowest_error, kwargs.get('lowest_error', 1.0))
        if kwargs.get('passed'):
            if kwargs.get('lowest_pass_error', self.lowest_pass_error) >= self.lowest_pass_error:
                return
            self.passed = True
            self.lowest_pass_error = kwargs['lowest_pass_error']
            self.validating_seq = kwargs.get('validating_seq')
            self.junction_results = kwargs.get('junction_results')
            self.source = kwargs.get('source')
            self.best_strand_match = kwargs.get('best_strand_match')

        elif not self.passed:
            if 'overlapping_spanning_reads_found' in kwargs:
                self.overlapping_spanning_reads_found = bool(kwargs['overlapping_spanning_reads_found'])
            if 'candidate_read_found' in kwargs:
                self.candidate_read_found = kwargs['candidate_read_found']
            if kwargs.get('junction_results'):
                self.junction_results = kwargs['junction_results']

        # Status / reason diagnostic updates based on new fields
        if self.source:
            if self.reference_check_match:
                self.status = SubseqStatus.INCONCLUSIVE
                self.reason = SubseqReason.REFERENCE_MATCH
            else:
                self.status = SubseqStatus.PASS
                self.reason = SubseqReason.PASS

        # Read based validation failure diagnostics
        elif not self.overlapping_spanning_reads_found:
            self.status = SubseqStatus.INCONCLUSIVE
            self.reason = SubseqReason.INCONCLUSIVE
        elif not self.candidate_read_found:
            self.status = SubseqStatus.FAIL
            self.reason = SubseqReason.OVER_ERROR_THRESHOLD

        # Assembly based validation failure diagnostics
        elif self.junction_results:
            self.status = SubseqStatus.FAIL
            self.reason = SubseqReason.JUNCTION_FAILED
        else:
            self.status = SubseqStatus.FAIL
            self.reason = SubseqReason.OTHER

    def apply_ambiguity(self, ambiguity: Optional['AmbiguityResult']) -> None:
        if ambiguity is not None and ambiguity.reference_match:
            self.reference_check_match = True
            self.reference_check_err = ambiguity.reference_err
            self.reference_check_strand = ambiguity.reference_strand

    def jsonify(self) -> Dict:
        return {
            'chrom': self.query.chrom,
            'ref_start': self.query.ref_start,
            'status': self.status,
            'reason': self.reason,
            'source': self.source,
            'strand': self.decisive_strand,
            'junctions': [j.jsonify() for j in self.junction_results],
            'validating_sequence': self.validating_seq,
        }

    @property
    def decisive_strand(self) -> Optional[int]:
        # Strand of the decisive alignment where meaningful: the reference match for a
        # reference_match, else the assembly match. None for read/edlib confirmations
        # (reads are randomly oriented; strand carries no signal there).
        return self.reference_check_strand if self.reason == SubseqReason.REFERENCE_MATCH else self.best_strand_match

     
@dataclass
class AmbiguityResult:
    reference_match: bool = False
    reference_err: float = 1.0
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
        # When True, assembly validation only considers forward-strand alignments
        # (reverse-strand hits are filtered out before the error/junction checks).
        self.assembly_forward_match_only = False

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

    def _assert_sv_records_contiguous_intervals(self, records: List[VariantRecord]) -> None:
        """
        Raises ValueError unless records of an sv satisfy:
        - all [start, stop) intervals are non overlapping and contiguous
        - all targets are either outside of [min(starts), max(stops)), or land on an existing start or stop
        """
        # sort by (start, stop), check for contiguity
        intervals = sorted(((*get_start_stop(rec), rec) for rec in records), key=lambda t: t[:2])
        for (prev_start, prev_stop, prev_rec), (start, stop, rec) in zip(intervals, intervals[1:]):
            if start != prev_stop:
                interval_error_type = 'overlap' if start < prev_stop else 'gap'
                raise ValueError(f'{rec.chrom}: record {rec.id} [{start},{stop}) is not contiguous with '
                                 f'the preceding record {prev_rec.id} [{prev_start},{prev_stop}) '
                                 f'({interval_error_type} of {abs(prev_stop - start)} bp)')

        # find min/max start, stops, check target pos
        starts = {start for start, _, _ in intervals}
        stops = {stop for _, stop, _ in intervals}
        min_start = intervals[0][0]
        max_stop = max(stop for _, stop, _ in intervals)
        for rec in records:
            target = rec.info.get('TARGET')
            if target is None:
                continue
            if min_start <= target < max_stop and target not in starts and target not in stops:
                raise ValueError(f'{rec.chrom}: record {rec.id} TARGET={target} falls inside the SV\'s '
                                 f'span [{min_start},{max_stop}) but does not land on an existing interval boundary')

    def _assert_sv_records_non_interchromosomal(self, records: List[VariantRecord]) -> None:
        """
        Raises ValueError if any record's TARGET_CHROM differs from the SV's own chromosome.
        """
        chrom = records[0].chrom
        for rec in records:
            target_chrom = rec.info.get('TARGET_CHROM', chrom)
            if target_chrom != chrom:
                raise ValueError(f'record {rec.id} on {chrom} has an interchromosomal TARGET_CHROM={target_chrom}')

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
                            'segments': [seg.jsonify() for seg in result.get('segments', [])],
                        }) + '\n')

                    total_calls[sv_type] += 1
                    overall_count += 1

                    if plot_first_n and total_calls[sv_type] <= plot_first_n:
                        query_infos = result.get('segments', [])
                        if query_infos:
                            logger.info(f'Producing dot plots for SV {svid} ({sv_type}, {len(query_infos)} part(s))')
                            for part, query_info in enumerate(query_infos):
                                plot_query_dot_plots(query_info, svid, sv_type, plot_out_dir, part=part)

                    if outcome == Outcome.HIT:
                        correct_calls[sv_type] += 1
                        overall_correct += 1
                        src = result.get('validation_source') or ValidationSource.ASSEMBLY
                        if src == ValidationSource.EDLIB:
                            edlib_rescues += 1
                        source_counts[src] += 1
                        if src == ValidationSource.READS:
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

            self._assert_sv_records_non_interchromosomal(records)
            self._assert_sv_records_contiguous_intervals(records)

            eval_reads = self.eval_mode in ('reads', 'both')
            eval_assembly = self.eval_mode in ('assembly', 'both')
            
            recon_sequences: List[QueryReconSubsequence] = simulate_subsequences(records, buffer, self.ref)
            score_records: List[QueryInfo] = []  # one per reconstructed subsequence of an SV


            candidate_reads_by_chrom = None
            if eval_reads:
                candidate_reads_by_chrom: Dict[str, list] = self.bam_reader.candidate_reads_from_records(records) if eval_reads else None # TODO: 

            for query in recon_sequences:
                # Window size 
                junction_validation_radius = int(clamp(round(self.junction_window_factor * len(query.sequence)), self.junction_window_min, self.junction_window_max))

                score_record = QueryInfo(query=query)
                if query.chrom in self.ref:
                    score_record.ref_sequence = self.ref[query.chrom][query.ref_start:query.ref_end].decode('ascii')

                if eval_reads:
                    read_based_score = self._score_read_based_eval(query, junction_validation_radius, candidate_reads_by_chrom[query.chrom], buffer)
                    score_record.update(**read_based_score)

                if not score_record.passed and eval_assembly:
                    assembly_based_score = self._score_assembly_based_eval(query, junction_validation_radius, location_tolerance, match_error_threshold)
                    score_record.update(**assembly_based_score)

                    if not score_record.passed and len(query.sequence) < MIN_EDLIB_QUERY:
                        edlib_based_score = self._score_edlib_based_eval(query, junction_validation_radius, match_error_threshold, buffer)
                        score_record.update(**edlib_based_score)

                if score_record.passed and self.check_reference:
                    score_record.apply_ambiguity(self._check_ambiguity(query, junction_validation_radius, match_error_threshold, buffer))
                
                score_records.append(score_record)

            coords = records[0].pos, records[0].stop

            # SV-level roll-up:
            #   hit          -> every subsequence passed
            #   miss         -> at least one subsequence was contradicted (spanning
            #                   reads / assembly aligned but didn't match)
            #   inconclusive -> no contradiction, but at least one subsequence could
            #                   not be tested (reads mode: no read spans the full allele)
            subseq_status = [r.status for r in score_records]
            if len(subseq_status) > 0 and all(s == SubseqStatus.PASS for s in subseq_status):
                outcome = Outcome.HIT
            elif SubseqStatus.FAIL in subseq_status:
                outcome = Outcome.MISS
            elif SubseqStatus.INCONCLUSIVE in subseq_status:
                outcome = Outcome.INCONCLUSIVE
            else:
                outcome = Outcome.MISS

            # Per-SV provenance: the highest tier any passing subsequence needed.
            # 'reads' implies the SV would have been an assembly-miss without reads.
            if outcome == Outcome.HIT:
                srcs = set(r.source for r in score_records if r.source)
                if ValidationSource.READS in srcs:
                    validation_source = ValidationSource.READS
                elif ValidationSource.EDLIB in srcs:
                    validation_source = ValidationSource.EDLIB
                else:
                    validation_source = ValidationSource.ASSEMBLY
            else:
                validation_source = None

            # Coarse evaluation tier for two-tier reporting:
            #   'read'     (Tier 1) -- a single actual read contained the allele. The
            #              strongest, most concrete evidence; independent of any
            #              assembly's quality.
            #   'assembly' (Tier 2) -- confirmed only against the reconstructed
            #              assembly (mappy or the edlib fallback). Supporting evidence
            #              for events too long for any single read to span.
            if validation_source == ValidationSource.READS:
                tier = Tier.READ
            elif validation_source in (ValidationSource.ASSEMBLY, ValidationSource.EDLIB):
                tier = Tier.ASSEMBLY
            else:
                tier = Tier.MISS

            return {
                'svid': svid,
                'sv_type': sv_type,
                'is_correct': outcome == Outcome.HIT,
                'outcome': outcome,   # 'hit' | 'miss' | 'inconclusive'
                'validation_source': validation_source,
                'tier': tier,   # 'read' (Tier 1) | 'assembly' (Tier 2) | None
                'coords': coords,
                'queries': [r.query for r in score_records],
                'match_scores': [r.lowest_pass_error for r in score_records],
                'best_scores': [r.lowest_error for r in score_records],
                'validating_sequences': [r.validating_seq for r in score_records],
                'recon_sequences': recon_sequences,
                'subseq_diagnostics': [r.reason.value for r in score_records],
                'segments': score_records,
            }
        
        except Exception as e:
            error_msg = f'Worker failed for SVID: {records[0].info.get("SVID", "Unknown")}. Error: {e}\n{traceback.format_exc()}'
            logger.error(error_msg)
            raise e
        
    def _score_read_based_eval(self, query: QueryReconSubsequence, junction_validation_radius, candidate_reads) -> Tuple[Dict, bool]:
        '''
        Read based validation,
        '''
        lowest_pass_error = 1.0
        lowest_error = 1.0
        passed = False
        validating_seq = None
        best_junction_validation_results = []
        candidate_read_found = False

        reads_spanning = [r for r in candidate_reads if len(r) >= len(query.sequence)]

        if reads_spanning:
            read_res: ReadEdlibResult = run_read_edlib(query.sequence, reads_spanning, self.read_error_threshold, self.min_read_support)
            candidate_read_found = read_res.was_read_found()
            if read_res.was_read_found():
                lowest_error = min(lowest_error, read_res.error)
                # TODO: check if this check is necessary
                if read_res.error <= self.read_error_threshold:
                    cigartuples = edlib_to_cigartuples(read_res.cigar)
                    junctions_validation_results: List[SeqSegmentValidationResult] = \
                                                                        validate_segments_from_cigar(cigartuples, query.segments,
                                                                                                radius=junction_validation_radius,
                                                                                                error_threshold=self.read_error_threshold)
                    if all(j.passed for j in junctions_validation_results):
                        passed = True
                        if read_res.error < lowest_pass_error:
                            lowest_pass_error = read_res.error
                            best_junction_validation_results = junctions_validation_results
                            validating_seq = read_res.matched_read_sequence
                    if not passed and read_res.error <= lowest_error:
                        best_junction_validation_results = junctions_validation_results

        return {
            'source': ValidationSource.READS,
            'lowest_pass_error': lowest_pass_error,
            'lowest_error': lowest_error,
            'passed': passed,
            'validating_seq': validating_seq,
            'junction_results': best_junction_validation_results,
            # read validation specific fields
            'overlapping_spanning_reads_found': reads_spanning,
            'candidate_read_found': candidate_read_found,
        }

    def _score_assembly_based_eval(self, query: QueryReconSubsequence, junction_validation_radius, location_tolerance, match_error_threshold):
        lowest_pass_error = 1.0
        lowest_error = 1.0
        passed = False
        validating_seq = None
        best_junction_validation_results = []
        # Assembly match specific fields
        best_strand_match = None

        for aligner in self.aligners[query.chrom]:
            all_alignments: List[mappy.Alignment] = list(aligner.map(query.sequence))
            alignments = [a for a in all_alignments if abs(a.r_st - query.ref_start) <= location_tolerance]
            if self.assembly_forward_match_only:
                alignments = [a for a in alignments if a.strand == 1]

            for alignment in alignments:
                match_err = check_match(alignment, query.sequence)
                lowest_error = min(lowest_error, match_err)
                if match_err <= match_error_threshold:
                    junctions_validation_results: List[SeqSegmentValidationResult] = \
                                                        validate_segments_from_cigar(alignment.cigar, query.segments,
                                                        radius=junction_validation_radius,
                                                        q_st=alignment.q_st, q_en=alignment.q_en,
                                                        query_len=len(query.sequence), strand=alignment.strand,
                                                        error_threshold=match_error_threshold)
                    if all(j.passed for j in junctions_validation_results):
                        passed = True
                        if match_err < lowest_pass_error:
                            lowest_pass_error = match_err
                            best_strand_match = alignment.strand   # record the best passing match's strand
                            best_junction_validation_results = junctions_validation_results
                            matched = aligner.seq(alignment.ctg, alignment.r_st, alignment.r_en)
                            validating_seq = reverse_complement(matched) if alignment.strand == -1 else matched
                    if not passed and match_err <= lowest_error:
                        best_junction_validation_results = junctions_validation_results

        return {
            'source': ValidationSource.ASSEMBLY,
            'lowest_pass_error': lowest_pass_error,
            'lowest_error': lowest_error,
            'passed': passed,
            'validating_seq': validating_seq,
            'junction_results': best_junction_validation_results,
            # assembly match specific fields
            'best_strand_match': best_strand_match,
        }

    def _score_edlib_based_eval(self, query: QueryReconSubsequence, junction_validation_radius, match_error_threshold, buffer):
        lowest_pass_error = 1.0
        lowest_error = 1.0
        passed = False
        validating_seq = None
        best_junction_validation_results = []

        for samp_bytes in self.sample:
            edlib_res = run_edlib_fallback(
                query.sequence,
                query.chrom,
                query.ref_start,
                buffer,
                EDLIB_FALLBACK_MAX_TOLERANCE,
                match_error_threshold,
                samp_bytes
            )
            if edlib_res:
                lowest_error = min(lowest_error, edlib_res.error)
            if edlib_res and edlib_res.error <= match_error_threshold:
                cigartuples = edlib_to_cigartuples(edlib_res.cigar)
                junctions_validation_results: List[SeqSegmentValidationResult] = \
                    validate_segments_from_cigar(cigartuples, query.segments,
                                                    radius=junction_validation_radius,
                                                    error_threshold=match_error_threshold)
                if all(j.passed for j in junctions_validation_results):
                    passed = True
                    if edlib_res.error < lowest_pass_error:
                        lowest_pass_error = edlib_res.error
                        validating_seq = edlib_res.matched_target_sequence
                        best_junction_validation_results = junctions_validation_results

                if not passed and edlib_res.error <= lowest_error:
                    best_junction_validation_results = junctions_validation_results

        return {
            'source': ValidationSource.EDLIB,
            'lowest_pass_error': lowest_pass_error,
            'lowest_error': lowest_error,
            'passed': passed,
            'validating_seq': validating_seq,
            'junction_results': best_junction_validation_results,
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
            ref_res = run_edlib_fallback(query.sequence, query.chrom, query.ref_start, int(buffer),
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