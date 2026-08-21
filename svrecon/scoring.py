"""Scoring engine: per-query/per-SV validation state and the CallsetScorer orchestrator."""
import json
import logging
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from pysam import VariantRecord
from tqdm import tqdm

from svrecon.constants import *
from svrecon.plot import plot_sv_validation
from svrecon.reconstruct import Query, simulate_subsequences
from svrecon.scorers.assembly import AssemblyScorer
from svrecon.scorers.base import Scorer, QueryValidationInput
from svrecon.scorers.edlib import EdlibScorer
from svrecon.scorers.reads import ReadScorer
from svrecon.utils import get_start_stop, group_variants_by_id, load_fasta_to_bytes

logger = logging.getLogger(__name__)


@dataclass
class QueryValidation:
    """Every scorer's result for one reconstructed subsequence, in the order they ran.
    Scoring stops at the first pass, so the last entry is the decisive one."""
    query: Query
    validations: List[QueryValidationInput] = field(default_factory=list)
    # Reference-ambiguity results, kept apart: a pass there means the allele also matches
    # the unmodified reference, making the query inconclusive rather than validated.
    ambiguity_validations: List[QueryValidationInput] = field(default_factory=list)

    def update_validation(self, result: QueryValidationInput) -> None:
        self.validations.append(result)

    def update_ambiguity(self, result: QueryValidationInput) -> None:
        self.ambiguity_validations.append(result)

    def jsonify(self) -> Dict:
        return {
            'chrom': self.query.chrom,
            'ref_start': self.query.ref_start,
            'status': self.status,
            'reason': self.reason,
            'source': self.source,
            'strand': self.decisive_strand,
            'validations': [v.jsonify() for v in self.validations],
            'ambiguity_validations': [v.jsonify() for v in self.ambiguity_validations],
        }

    @property
    def latest(self) -> Optional[QueryValidationInput]:
        return self.validations[-1] if self.validations else None

    @property
    def passed(self) -> bool:
        return bool(self.latest and self.latest.passed)

    @property
    def reference_matches(self) -> List[QueryValidationInput]:
        return [v for v in self.ambiguity_validations if v.passed]

    @property
    def status(self) -> SubseqStatus:
        if self.reference_matches:
            return SubseqStatus.INCONCLUSIVE
        return self.latest.status if self.latest else SubseqStatus.SKIPPED

    @property
    def reason(self) -> SubseqReason:
        if self.reference_matches:
            return SubseqReason.REFERENCE_MATCH
        return self.latest.reason if self.latest else SubseqReason.SKIPPED

    @property
    def source(self) -> Optional[ValidationSource]:
        return self.latest.source if self.passed else None

    @property
    def validating_seq(self) -> Optional[str]:
        return self.latest.validating_seq if self.passed else None

    @property
    def lowest_pass_error(self) -> float:
        return self.latest.lowest_pass_error if self.passed else 1.0

    @property
    def decisive_strand(self) -> Optional[int]:
        # Strand of the alignment that decided the query: the reference match if the allele
        # also matched the reference, else the validating one. None for reads (randomly
        # oriented, so strand carries no signal).
        if not self.passed:
            return None
        return (self.reference_matches or self.validations)[-1].best_strand_match


class SVValidation:
    """SV-level roll-up of one call's per-subsequence validations."""

    def __init__(self, svid: str, svtype: str, query_validations: List[QueryValidation]):
        self.svid = svid
        self.svtype = svtype
        self.query_validations = query_validations

        # SV-level roll-up:
        #   skipped      -> no validation scorers configured: reconstructed but never
        #                   validated at all -- distinct from "inconclusive" (which
        #                   means validation was attempted but couldn't reach a verdict)
        #   hit          -> every subsequence passed
        #   miss         -> at least one subsequence was contradicted (spanning
        #                   reads / assembly aligned but didn't match)
        #   inconclusive -> no contradiction, but at least one subsequence could
        #                   not be tested (reads mode: no read spans the full allele)
        subseq_status = [qv.status for qv in query_validations]
        if len(subseq_status) > 0 and all(s == SubseqStatus.SKIPPED for s in subseq_status):
            self.outcome = Outcome.SKIPPED
        elif len(subseq_status) > 0 and all(s == SubseqStatus.PASS for s in subseq_status):
            self.outcome = Outcome.HIT
        elif SubseqStatus.FAIL in subseq_status:
            self.outcome = Outcome.MISS
        elif SubseqStatus.INCONCLUSIVE in subseq_status:
            self.outcome = Outcome.INCONCLUSIVE
        else:
            self.outcome = Outcome.MISS

        # Per-SV provenance: the highest tier any passing subsequence needed.
        # 'reads' implies the SV would have been an assembly-miss without reads.
        self.validation_source = None
        if self.outcome == Outcome.HIT:
            srcs = {qv.source for qv in query_validations if qv.source}
            if ValidationSource.READS in srcs:
                self.validation_source = ValidationSource.READS
            elif ValidationSource.EDLIB in srcs:
                self.validation_source = ValidationSource.EDLIB
            else:
                self.validation_source = ValidationSource.ASSEMBLY

        # Coarse evaluation tier for two-tier reporting:
        #   'read'     (Tier 1) -- a single actual read contained the allele. The
        #              strongest, most concrete evidence; independent of any
        #              assembly's quality.
        #   'assembly' (Tier 2) -- confirmed only against the reconstructed
        #              assembly (mappy or the edlib fallback). Supporting evidence
        #              for events too long for any single read to span.
        if self.validation_source == ValidationSource.READS:
            self.tier = Tier.READ
        elif self.validation_source in (ValidationSource.ASSEMBLY, ValidationSource.EDLIB):
            self.tier = Tier.ASSEMBLY
        else:
            self.tier = Tier.MISS

    def get_summary(self) -> Dict:
        """Report record for the --report json sidecar."""
        return {
            'svid': self.svid,
            'svtype': self.svtype,
            'outcome': self.outcome,
            'tier': self.tier,
            'query_validations': [qv.jsonify() for qv in self.query_validations],
        }

    @property
    def chroms(self) -> List[str]:
        """Chromosome of each reconstructed subsequence (an SV can span multiple)."""
        return [qv.query.chrom for qv in self.query_validations]

    @property
    def match_scores(self) -> List[float]:
        return [qv.lowest_pass_error for qv in self.query_validations]

    @property
    def subseq_reasons(self) -> List[str]:
        return [qv.reason.value for qv in self.query_validations]


class CallsetScorer(object):
    def __init__(self, config):
        # Paths & plotting
        self.report_path = config.report_path
        self.plot_first_n = config.plot_first_n
        self.plot_out_dir = config.img_dir
        self.plot_aspect = config.plot_aspect

        # Parse and group records by svid
        self.variants = group_variants_by_id(config.calls, config.gap_file)
        logger.info(f'Found {len(self.variants)} calls from callset {config.calls}')

        # Identify chromsomes of interest in callset
        chroms = set()
        for records in self.variants.values():
            for record in records:
                chroms.add(record.chrom)
                if 'TARGET_CHROM' in record.info:
                    chroms.add(record.info['TARGET_CHROM'])
        logger.info(f'Found {len(chroms)} referenced chromosomes in callset: {chroms}')

        # Shared parameters; reference bytes for reconstruction
        self.buffer = config.buffer
        self.auto_buffer_min = config.auto_buffer_min
        self.auto_buffer_fraction = config.auto_buffer_fraction
        self.n_threads = config.n_threads
        self.verbose = config.verbose
        self.chrom_to_ref = load_fasta_to_bytes(config.reference, chroms)
        logger.info(f'Loaded {len(self.chrom_to_ref)} reference chromosomes from {config.reference}')

        # Populate Validation Scorers, sorted by tier
        self.validation_scorer: list[Scorer] = []
        self.ambiguity_scorers: list[Scorer] = []

        if config.bam:
            self.validation_scorer.append(ReadScorer(config, chroms))  # TODO: ReadScorer.score_query pending
            logger.info("Initialized read scoring validation")

        if config.sample:
            self.validation_scorer.extend([AssemblyScorer(config, chroms, config.sample, config.assembly_forward_match_only),
                                           EdlibScorer(config, chroms, config.sample)])
            logger.info("Initialized assembly scoring validation")

        # Populate Ambiguity Scorers: the same machinery aimed at the reference -- a pass there
        # means the allele also matches the unmodified reference (not specific to the SV).
        if config.check_reference:
            self.ambiguity_scorers.extend([AssemblyScorer(config, chroms, config.reference, forward_match_only=False),
                                           EdlibScorer(config, chroms, config.reference)])
            logger.info("Initialized reference ambiguity scoring")

    def _resolve_buffer(self, records: List[VariantRecord]) -> int:
        """Resolves the 'auto' sentinel to max(auto_buffer_min, auto_buffer_fraction * the SV's
        longest segment); a fixed buffer passes through unchanged."""
        if self.buffer != 'auto':
            return int(self.buffer)
        segment_lengths = [stop - start for start, stop in {get_start_stop(rec) for rec in records}]
        if not segment_lengths:
            return self.auto_buffer_min
        return max(self.auto_buffer_min, int(self.auto_buffer_fraction * max(segment_lengths)))

    def _assert_sv_records_contiguous_intervals(self, records: List[VariantRecord]) -> None:
        """
        Asserts the following conditions:
        - all [start, stop) intervals are non overlapping and contiguous
        - all targets are either outside of [min(starts), max(stops)], or land on an existing endpoint
        """
        # Dedupe identical (start, stop) spans -- multiple records can share one source span
        # for different downstream operations -- before sorting for contiguity.
        by_interval: Dict[Tuple[int, int], VariantRecord] = {}
        for rec in records:
            by_interval.setdefault(get_start_stop(rec), rec)
        intervals = sorted((start, stop, rec) for (start, stop), rec in by_interval.items())
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
            if min_start <= target <= max_stop and target not in starts and target not in stops:
                raise ValueError(f'{rec.chrom}: record {rec.id} TARGET={target} falls inside the SV\'s '
                                 f'span [{min_start},{max_stop}) but does not land on an existing interval boundary')

    def _assert_sv_records_non_interchromosomal(self, records: List[VariantRecord]) -> None:
        """
        Raises ValueError if records contains more than one unique chromosome
        """
        chrom = records[0].chrom
        for rec in records:
            target_chrom = rec.info.get('TARGET_CHROM', chrom)
            if target_chrom != chrom:
                raise ValueError(f'record {rec.id} on {chrom} has an interchromosomal TARGET_CHROM={target_chrom}')

    def score_all(self):
        total_calls = Counter()
        correct_calls = Counter()
        inconclusive_calls = Counter()  # reads mode: no read spans the full resulting allele
        skipped_calls = Counter()   # neither --sample nor --bam: reconstructed but never validated
        assembly_hits = Counter()   # hits validated by the assembly (mappy or edlib)
        read_hits = Counter()       # hits validated only by reads
        source_counts = Counter()   # overall hit provenance: assembly | edlib | reads
        overall_count = 0
        overall_correct = 0
        overall_inconclusive = 0
        overall_skipped = 0
        edlib_rescues = 0
        precision = {}

        # Optional per-SV JSON sidecar (one record per line). Additive: it does not affect the
        # log output or the score table. Written from this single consumer thread, so no lock.
        report_fh = open(self.report_path, 'w') if self.report_path else None

        with ThreadPoolExecutor(max_workers=self.n_threads) as executor:
            futures = {executor.submit(self.score_sv, self.variants[svid]) for svid in self.variants.keys()}

            pbar = tqdm(as_completed(futures), total=len(self.variants), desc=f'Scoring SVs', smoothing=0)

            for future in pbar:
                sv_validation: SVValidation = future.result()
                svid, sv_type, outcome = sv_validation.svid, sv_validation.svtype, sv_validation.outcome
                line = f'{svid}\t{sv_type}\t{sv_validation.chroms}\t{outcome.value}\t{sv_validation.match_scores}'
                if outcome == Outcome.HIT:
                    # Tier 1 (read) vs Tier 2 (assembly), plus the detailed source.
                    line += f'\t{sv_validation.tier}\t{sv_validation.validation_source}'
                elif outcome != Outcome.SKIPPED:  # miss or inconclusive -> show why
                    line += f'\t{sv_validation.subseq_reasons}'
                if self.verbose:
                    line += f'\t{json.dumps(sv_validation.get_summary())}'
                logger.info(line)

                if report_fh is not None:
                    report_fh.write(json.dumps(sv_validation.get_summary()) + '\n')

                total_calls[sv_type] += 1
                overall_count += 1

                if self.plot_first_n and total_calls[sv_type] <= self.plot_first_n:
                    if sv_validation.query_validations:
                        logger.info(f'Producing dot plots for SV {svid} ({sv_type}, '
                                    f'{len(sv_validation.query_validations)} part(s))')
                        outcome_dir = 'validated' if outcome == Outcome.HIT else 'unvalidated'
                        sv_plot_dir = Path(self.plot_out_dir) / outcome_dir / sv_type / svid
                        sv_plot_dir.mkdir(parents=True, exist_ok=True)
                        plot_sv_validation(sv_validation, str(sv_plot_dir), aspect=self.plot_aspect)

                if outcome == Outcome.HIT:
                    correct_calls[sv_type] += 1
                    overall_correct += 1
                    src = sv_validation.validation_source or ValidationSource.ASSEMBLY
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
                elif outcome == Outcome.SKIPPED:
                    skipped_calls[sv_type] += 1
                    overall_skipped += 1
                # else: miss -> counts toward total but not correct/inconclusive/skipped

                # precision is over CONCLUSIVE calls only (hits + misses; inconclusive and
                # skipped -- never validated at all -- are excluded from the denominator)
                conclusive = overall_count - overall_inconclusive - overall_skipped
                if conclusive:
                    desc = (f'Scoring SVs. Precision {overall_correct / conclusive:.2f} '
                           f'({overall_correct}/{conclusive}); inconclusive {overall_inconclusive}')
                    if overall_skipped:
                        desc += f'; skipped {overall_skipped}'
                else:
                    desc = f'Scoring SVs. {overall_skipped} skipped so far, 0 conclusive'
                pbar.set_description(desc)
                denom = total_calls[sv_type] - inconclusive_calls[sv_type] - skipped_calls[sv_type]
                precision[sv_type] = correct_calls[sv_type] / denom if denom else 0

        if report_fh is not None:
            report_fh.close()
            logger.info(f'Wrote per-SV eval report: {self.report_path}')

        logger.info(f'Total SVs rescued by edlib fallback: {edlib_rescues}')
        overall_miss = overall_count - overall_correct - overall_inconclusive - overall_skipped
        logger.info(
            f"Outcomes: hit={overall_correct} miss={overall_miss} inconclusive={overall_inconclusive} "
            f"skipped={overall_skipped} "
            f"(inconclusive = no read spans the full resulting allele; skipped = neither --sample "
            f"nor --bam given; both excluded from precision)")
        logger.info(
            f"Hits by validation source: assembly={source_counts['assembly']} "
            f"edlib={source_counts['edlib']} reads={source_counts['reads']}")
        logger.info(
            f"Hits by tier: read(Tier1)={source_counts['reads']} "
            f"assembly(Tier2)={source_counts['assembly'] + source_counts['edlib']} "
            f"(Tier1 = a single real read contained the allele; Tier2 = confirmed only against the assembly)")

        total_all = sum(total_calls.values())
        inc_all = sum(inconclusive_calls.values())
        skip_all = sum(skipped_calls.values())
        conclusive_all = total_all - inc_all - skip_all
        precision['ALL'] = sum(correct_calls.values()) / conclusive_all if conclusive_all else 0
        correct_calls['ALL'] = sum(correct_calls.values())
        total_calls['ALL'] = total_all
        inconclusive_calls['ALL'] = inc_all
        skipped_calls['ALL'] = skip_all
        assembly_hits['ALL'] = sum(assembly_hits.values())
        read_hits['ALL'] = sum(read_hits.values())

        return precision, correct_calls, total_calls, inconclusive_calls, skipped_calls, assembly_hits, read_hits

    def score_sv(self, records: List[VariantRecord]) -> SVValidation:
        """Score one SV: reconstruct its alt allele(s), validate each subsequence with the
        configured scorers, and return the SVValidation."""
        svid = records[0].info['SVID']
        sv_type = records[0].info['SVTYPE']

        try:
            self._assert_sv_records_non_interchromosomal(records)
            self._assert_sv_records_contiguous_intervals(records)

            sv_buffer = self._resolve_buffer(records)
            query_validations: List[QueryValidation] = []  # one per reconstructed subsequence of an SV

            for query in simulate_subsequences(records, sv_buffer, self.chrom_to_ref):
                query_validation = QueryValidation(query=query)

                for scorer in self.validation_scorer:
                    query_validation.update_validation(scorer.score_query(query))
                    if query_validation.passed: break

                if query_validation.passed:
                    for scorer in self.ambiguity_scorers:
                        query_validation.update_ambiguity(scorer.score_query(query))
                        if query_validation.reference_matches: break

                query_validations.append(query_validation)

            return SVValidation(svid, sv_type, query_validations)
        except Exception as e:
            e.add_note(f'while scoring SV {svid} ({sv_type})')
            raise
