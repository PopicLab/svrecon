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
from svrecon.io import write_annotated_vcf
from svrecon.plot import plot_sv_validation
from svrecon.reconstruct import Query, construct_queries
from svrecon.scorers.assembly import AssemblyScorer
from svrecon.scorers.base import Scorer, QueryValidation
from svrecon.scorers.edlib import EdlibScorer
from svrecon.scorers.reads import ReadScorer
from svrecon.utils import get_start_stop, load_grouped_variants_from_vcf, load_fasta_to_bytes

logger = logging.getLogger(__name__)


@dataclass
class QueryValidationResult:
    """Validation results per reconstructed query summarized over cigar checks"""
    query: Query
    validations: List[QueryValidation] = field(default_factory=list)
    ambiguity_validations: List[QueryValidation] = field(default_factory=list)

    def update_validation(self, result: QueryValidation) -> None:
        self.validations.append(result)

    def update_ambiguity(self, result: QueryValidation) -> None:
        self.ambiguity_validations.append(result)

    def jsonify(self) -> Dict:
        return {
            'query': self.query.jsonify(),
            'status': self.status,
            'reason': self.reason,
            'aligned': self.aligned,
            'reference_ambiguous': bool(self.reference_matches),
            'source': self.source,
            'strand': self.decisive_strand,
            'validations': [v.jsonify() for v in self.validations],
            'ambiguity_validations': [v.jsonify() for v in self.ambiguity_validations],
        }

    @property
    def decisive(self) -> Optional[QueryValidation]:
        """The result that settles the query -- see QueryValidation.__gt__."""
        return max(self.validations) if self.validations else None

    @property
    def passed(self) -> bool:
        return bool(self.decisive and self.decisive.passed)

    @property
    def reference_matches(self) -> List[QueryValidation]:
        return [v for v in self.ambiguity_validations if v.passed]

    @property
    def status(self) -> QueryValidationStatus:
        # No decisive result means no scorer was configured -- untestable, like a reference match.
        if self.reference_matches or not self.decisive:
            return QueryValidationStatus.INCONCLUSIVE
        return self.decisive.status

    @property
    def reason(self) -> QueryValidationReason:
        if self.reference_matches:
            return QueryValidationReason.REFERENCE_AMBIGUOUS
        if not self.decisive:
            return QueryValidationReason.NOT_ATTEMPTED
        return self.decisive.reason

    @property
    def aligned(self) -> bool:
        return bool(self.decisive and self.decisive.aligned)

    @property
    def source(self) -> Optional[ValidationSource]:
        return self.decisive.source if self.decisive else None

    @property
    def best_matched_seq(self) -> Optional[str]:
        return self.decisive.best_matched_seq if self.decisive else None

    @property
    def lowest_pass_error(self) -> float:
        return self.decisive.lowest_pass_error if self.passed else 1.0

    @property
    def decisive_strand(self) -> Optional[int]:
        # Strand of the alignment that decided the query: the reference match if the allele
        # also matched the reference, else the validating one. None for reads (randomly
        # oriented, so strand carries no signal).
        if not self.passed:
            return None
        return self.reference_matches[-1].best_strand_match if self.reference_matches \
            else self.decisive.best_strand_match


class SVValidationResult:
    """SV-level roll-up of one call's per-query validations."""

    def __init__(self, svid: str, svtype: str, query_validations_results: List[QueryValidationResult], records: list[VariantRecord]):
        self.svid = svid
        self.svtype = svtype
        self.query_validation_results = query_validations_results
        self.records = records

        # SV-level roll-up:
        #   hit          -> every query passed
        #   miss         -> at least one query was contradicted (spanning reads /
        #                   assembly aligned but didn't match, or nothing aligned)
        #   inconclusive -> no contradiction, but at least one query could not be
        #                   tested: no read spans the full allele, the allele also
        #                   matches the reference, or no scorer was configured at all
        query_statuses = [qvr.status for qvr in query_validations_results]
        if query_statuses and all(s is QueryValidationStatus.PASS for s in query_statuses):
            self.outcome = Outcome.HIT
        elif any(s is QueryValidationStatus.FAIL for s in query_statuses):
            self.outcome = Outcome.MISS
        elif QueryValidationStatus.INCONCLUSIVE in query_statuses:
            self.outcome = Outcome.INCONCLUSIVE
        else:
            self.outcome = Outcome.MISS

        # Per-SV provenance: the highest tier any passing query needed.
        # 'reads' implies the SV would have been an assembly-miss without reads.
        self.validation_source = None
        if self.outcome == Outcome.HIT:
            srcs = {qv.source for qv in query_validations_results if qv.source}
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
            'query_validations': [qv.jsonify() for qv in self.query_validation_results],
        }

    @property
    def chroms(self) -> List[str]:
        """Chromosome of each reconstructed query (an SV can span multiple)."""
        return [qv.query.chrom for qv in self.query_validation_results]

    @property
    def match_scores(self) -> List[float]:
        return [qv.lowest_pass_error for qv in self.query_validation_results]

    @property
    def query_diagnostics(self) -> List[str]:
        """Per-query '<status>:<reason>' tags for the non-hit log lines."""
        return [f'{qv.status.value}:{qv.reason.value}' for qv in self.query_validation_results]


class CallsetScorer(object):
    def __init__(self, config):
        # Paths & plotting
        self.report_path = config.report_path
        self.annotated_vcf_path = config.annotated_vcf_path
        self.plot_first_n = config.plot_first_n
        self.plot_out_dir = config.img_dir
        self.plot_aspect = config.plot_aspect
        self.plot_substitute_bases = config.plot_substitute_bases
        self.plot_title_svid = config.plot_title_svid
        self.plot_title_location = config.plot_title_location
        self.plot_axis_length = config.plot_axis_length

        # Parse and group records by svid
        self.variants = load_grouped_variants_from_vcf(config.calls, config.gap_file)
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

    def _assert_sv_records_positive_lengths(self, records: List[VariantRecord]) -> None:
        """
        Asserts every [start, stop) interval is non empty
        """
        for rec in records:
            start, stop = get_start_stop(rec)
            if stop <= start:
                raise ValueError(f'{rec.chrom}: record {rec.id} interval [{start},{stop}) is empty; '
                                 f'stop must exceed start')

    def _assert_sv_records_non_overlapping(self, records: List[VariantRecord]) -> None:
        """
        Asserts all [start, stop) intervals are non overlapping
        """
        # Dedupe identical (start, stop) spans -- ensure non overlapping
        by_interval: Dict[Tuple[int, int], VariantRecord] = {}
        for rec in records:
            by_interval.setdefault(get_start_stop(rec), rec)
        intervals = sorted((start, stop, rec) for (start, stop), rec in by_interval.items())
        for (prev_start, prev_stop, prev_rec), (start, stop, rec) in zip(intervals, intervals[1:]):
            if start < prev_stop:
                raise ValueError(f'{rec.chrom}: record {rec.id} [{start},{stop}) overlaps '
                                 f'the preceding record {prev_rec.id} [{prev_start},{prev_stop}) '
                                 f'by {prev_stop - start} bp')

    def _assert_sv_records_no_split_segments(self, records: List[VariantRecord]) -> None:
        """
        Asserts no target lands in an interval's (start, stop) -- which would result in a split segment
        """
        for rec in records:
            if 'TARGET' not in rec.info:
                continue
            target = rec.info['TARGET']
            for source_rec in records:
                start, stop = get_start_stop(source_rec)
                if start < target < stop:
                    raise ValueError(f'{rec.chrom}: record {rec.id} TARGET={target} falls inside '
                                     f'record {source_rec.id}\'s interval [{start},{stop})')

    def _assert_sv_records_required_fields(self, records: List[VariantRecord]) -> None:
        """
        Asserts every record has an SVTYPE, and every record of a multi-record SV an OP_TYPE
        """
        required = ('SVTYPE',) if len(records) == 1 else ('SVTYPE', 'OP_TYPE')
        for rec in records:
            for key in required:
                if key not in rec.info:
                    raise ValueError(f'{rec.chrom}: record {rec.id} is missing INFO/{key}')
        

    def _assert_sv_records_non_interchromosomal(self, records: List[VariantRecord]) -> None:
        """
        Raises ValueError if records contains more than one unique chromosome
        """
        chrom = records[0].chrom
        for rec in records:
            target_chrom = rec.info['TARGET_CHROM'] if 'TARGET_CHROM' in rec.info else chrom
            if target_chrom != chrom:
                raise ValueError(f'record {rec.id} on {chrom} has an interchromosomal TARGET_CHROM={target_chrom}')

    def score_all(self):
        total_calls = Counter()
        correct_calls = Counter()
        inconclusive_calls = Counter()  # untestable: no spanning read, reference-ambiguous, or no scorer
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
        report_fh = open(self.report_path, 'w') if self.report_path else None
        svid_to_result: Dict[str, SVValidationResult] = {}

        with ThreadPoolExecutor(max_workers=self.n_threads) as executor:
            futures = {executor.submit(self.score_sv, svid, records)
                       for svid, records in self.variants.items()}

            pbar = tqdm(as_completed(futures), total=len(self.variants), desc=f'Scoring SVs', smoothing=0)

            for future in pbar:
                sv_validation: SVValidationResult = future.result()
                svid, sv_type, outcome = sv_validation.svid, sv_validation.svtype, sv_validation.outcome
                svid_to_result[svid] = sv_validation
                line = f'{svid}\t{sv_type}\t{sv_validation.chroms}\t{outcome.value}\t{sv_validation.match_scores}'
                if outcome == Outcome.HIT:
                    # Tier 1 (read) vs Tier 2 (assembly), plus the detailed source.
                    line += f'\t{sv_validation.tier}\t{sv_validation.validation_source}'
                else:  # miss or inconclusive -> show why
                    line += f'\t{sv_validation.query_diagnostics}'
                if self.verbose:
                    line += f'\t{json.dumps(sv_validation.get_summary(), default=str, indent=2)}'
                logger.info(line)

                if report_fh is not None:
                    report_fh.write(json.dumps(sv_validation.get_summary(), default=str) + '\n')

                total_calls[sv_type] += 1
                overall_count += 1

                if self.plot_first_n and total_calls[sv_type] <= self.plot_first_n:
                    if sv_validation.query_validation_results:
                        logger.info(f'Producing dot plots for SV {svid} ({sv_type}, '
                                    f'{len(sv_validation.query_validation_results)} part(s))')
                        outcome_dir = 'validated' if outcome == Outcome.HIT else 'unvalidated'
                        sv_plot_dir = Path(self.plot_out_dir) / outcome_dir / sv_type / svid
                        sv_plot_dir.mkdir(parents=True, exist_ok=True)
                        plot_sv_validation(sv_validation, str(sv_plot_dir), aspect=self.plot_aspect,
                                           substitute_bases=self.plot_substitute_bases,
                                           title_svid=self.plot_title_svid,
                                           title_location=self.plot_title_location,
                                           axis_length=self.plot_axis_length)

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
                # else: miss -> counts toward total but not correct/inconclusive

                # precision is over CONCLUSIVE calls only (hits + misses; inconclusive calls,
                # which no scorer could test, are excluded from the denominator)
                conclusive = overall_count - overall_inconclusive
                if conclusive:
                    desc = (f'Scoring SVs. Precision {overall_correct / conclusive:.2f} '
                           f'({overall_correct}/{conclusive}); inconclusive {overall_inconclusive}')
                else:
                    desc = f'Scoring SVs. {overall_inconclusive} inconclusive so far, 0 conclusive'
                pbar.set_description(desc)
                denom = total_calls[sv_type] - inconclusive_calls[sv_type]
                precision[sv_type] = correct_calls[sv_type] / denom if denom else 0

        sv_validations = [svid_to_result[svid] for svid in self.variants] # reload in original order by processed svid
        write_annotated_vcf(sv_validations, str(self.annotated_vcf_path))
        logger.info(f'Wrote annotated VCF: {self.annotated_vcf_path}')

        if report_fh is not None:
            report_fh.close()
            logger.info(f'Wrote per-SV eval report: {self.report_path}')

        logger.info(f'Total SVs rescued by edlib fallback: {edlib_rescues}')
        overall_miss = overall_count - overall_correct - overall_inconclusive
        logger.info(
            f"Outcomes: hit={overall_correct} miss={overall_miss} inconclusive={overall_inconclusive} "
            f"(inconclusive = no scorer could test the call -- no read spans the full resulting "
            f"allele, the allele also matches the reference, or neither --sample nor --bam was "
            f"given; excluded from precision)")
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
        precision['ALL'] = sum(correct_calls.values()) / conclusive_all if conclusive_all else 0
        correct_calls['ALL'] = sum(correct_calls.values())
        total_calls['ALL'] = total_all
        inconclusive_calls['ALL'] = inc_all
        assembly_hits['ALL'] = sum(assembly_hits.values())
        read_hits['ALL'] = sum(read_hits.values())

        return precision, correct_calls, total_calls, inconclusive_calls, assembly_hits, read_hits

    def score_sv(self, svid: str, records: List[VariantRecord]) -> SVValidationResult:
        """Score one SV: reconstruct its alt allele(s), validate each query with the
        configured scorers, and return the SVValidation. ``svid`` is the grouping key"""
        sv_type = records[0].info.get('SVTYPE')  # asserted present below, before any use

        try:
            self._assert_sv_records_required_fields(records)
            self._assert_sv_records_non_interchromosomal(records)
            self._assert_sv_records_positive_lengths(records)
            self._assert_sv_records_non_overlapping(records)
            self._assert_sv_records_no_split_segments(records)

            sv_buffer = self._resolve_buffer(records)
            query_validations: List[QueryValidationResult] = []  # one per reconstructed query of an SV

            for query in construct_queries(records, sv_buffer, self.chrom_to_ref):
                query_validation = QueryValidationResult(query=query)

                for scorer in self.validation_scorer:
                    query_validation.update_validation(scorer.score_query(query))
                    if query_validation.passed: break

                if query_validation.passed:
                    for scorer in self.ambiguity_scorers:
                        query_validation.update_ambiguity(scorer.score_query(query))
                        if query_validation.reference_matches: break

                query_validations.append(query_validation)

            return SVValidationResult(svid, sv_type, query_validations, records) 
        except Exception as e:
            raise RuntimeError(f'while scoring SV {svid} ({sv_type})') from e
