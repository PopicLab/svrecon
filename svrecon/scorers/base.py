"""Scorer classes: per-query validation against reads, a sample assembly (mappy), and edlib windows."""
import logging
from dataclasses import dataclass, field
from typing import Dict, List, Optional

import mappy

from svrecon.config import Config
from svrecon.constants import Outcome, SubseqReason, SubseqStatus, Tier, ValidationSource
from svrecon.reconstruct import Query
from svrecon.scorers.align import (build_chrom_aligners, check_match, run_edlib_fallback,
                                   validate_segments_from_cigar, SeqSegmentValidationResult)
from svrecon.scorers.reads import BamReader, run_read_edlib
from svrecon.scorers.utils import Cigar
from svrecon.util import load_fasta_to_bytes, reverse_complement

logger = logging.getLogger(__name__)


@dataclass
class QueryValidationInput:
    """One scorer's lightweight result for one query; QueryValidation folds these in."""
    source: ValidationSource
    passed: bool = False
    lowest_pass_error: float = 1.0
    lowest_error: float = 1.0
    validating_seq: Optional[str] = None
    segment_results: List[SeqSegmentValidationResult] = field(default_factory=list)
    best_strand_match: Optional[int] = None        # assembly only
    overlapping_spanning_reads_found: bool = True  # reads only; True = n/a for other tiers
    candidate_read_found: bool = True              # reads only


@dataclass
class QueryValidation:
    """Validation state of one reconstructed subsequence, accumulated across scorers'
    QueryValidationInput results."""
    query: Query
    source: Optional[ValidationSource] = None

    # Pass Diagnostics
    lowest_pass_error: float = 1.0
    lowest_error: float = 1.0  # lowest error seen overall, pass or fail
    passed: bool = False
    validating_seq: Optional[str] = None
    segment_results: List[SeqSegmentValidationResult] = field(default_factory=list)
    best_strand_match: Optional[int] = None  # assembly only

    # Failure diagnostics -- only meaningful once a tier call did NOT pass. Only the
    # reads tier emits these; default True ("not applicable/assume satisfied") so an
    # assembly/edlib-only call never spuriously trips the classification below.
    overlapping_spanning_reads_found: bool = True
    candidate_read_found: bool = True

    # Populated by apply_ambiguity() -- only meaningful when passed, and only checked
    # when check_reference is on.
    reference_check_match: bool = False
    reference_check_err: Optional[float] = None
    reference_check_strand: Optional[int] = None

    status: SubseqStatus = SubseqStatus.SKIPPED # by default, unscored
    reason: SubseqReason = SubseqReason.SKIPPED

    def update_validation(self, result: QueryValidationInput) -> None:
        """Fold one scorer's result in: a pass adopts the result's fields; failure
        diagnostics accumulate only while nothing has passed."""
        self.lowest_error = min(self.lowest_error, result.lowest_error)
        if result.passed:
            self.passed = True
            self.lowest_pass_error = result.lowest_pass_error
            self.validating_seq = result.validating_seq
            self.segment_results = result.segment_results
            self.source = result.source
            self.best_strand_match = result.best_strand_match

        elif not self.passed:
            if result.source == ValidationSource.READS:  # only the reads tier emits these flags
                self.overlapping_spanning_reads_found = result.overlapping_spanning_reads_found
                self.candidate_read_found = result.candidate_read_found
            if result.segment_results:
                self.segment_results = result.segment_results

        self._update_status()

    def _update_status(self) -> None:
        # Status / reason diagnostics based on the fields accumulated so far
        if self.source:
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
        elif self.segment_results:
            self.status = SubseqStatus.FAIL
            self.reason = SubseqReason.JUNCTION_FAILED
        else:
            self.status = SubseqStatus.FAIL
            self.reason = SubseqReason.OTHER

    def update_ambiguity(self, result: QueryValidationInput) -> None:
        """A pass against the plain reference is not specific to the SV -> inconclusive."""
        if result.passed:
            self.reference_check_match = True
            self.reference_check_err = result.lowest_pass_error
            self.reference_check_strand = result.best_strand_match or 1  # edlib aligns forward only
            self.status = SubseqStatus.INCONCLUSIVE
            self.reason = SubseqReason.REFERENCE_MATCH

    def jsonify(self) -> Dict:
        return {
            'chrom': self.query.chrom,
            'ref_start': self.query.ref_start,
            'status': self.status,
            'reason': self.reason,
            'source': self.source,
            'strand': self.decisive_strand,
            'segments': [s.jsonify() for s in self.segment_results],
            'validating_sequence': self.validating_seq,
        }

    @property
    def decisive_strand(self) -> Optional[int]:
        # Strand of the decisive alignment where meaningful: the reference match for a
        # reference_match, else the assembly match. None for read/edlib confirmations
        # (reads are randomly oriented; strand carries no signal there).
        return self.reference_check_strand if self.reason == SubseqReason.REFERENCE_MATCH else self.best_strand_match


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
            'segments': [qv.jsonify() for qv in self.query_validations],
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


class Scorer:
    """Base for per-query validators; subclasses implement score_query(query) -> QueryValidationInput."""

    def __init__(self, config: Config, chroms: set[str]):
        # general shared params
        self.match_error_threshold = config.match_error_threshold # for both reference and assembly alignment
        self.location_tolerance = config.location_tolerance
        self.chroms = chroms

    def score_query(self, query: Query) -> QueryValidationInput:
        raise NotImplementedError("Implement in subclass")


class ReadScorer(Scorer):
    """Validates a query against real reads: does any single read contain the allele?"""

    def __init__(self, config, chroms):
        super().__init__(config, chroms)
        self.read_error_threshold = config.read_error_threshold
        self.min_read_support = config.min_read_support
        self.max_reads_per_site = config.max_reads_per_site
        self.bam_reader = BamReader(config.bam)

    def score_query(self, query: Query) -> QueryValidationInput:
        # TODO: logic so far, to finish and re-enable:
        lowest_pass_error = 1.0
        lowest_error = 1.0
        passed = False
        validating_seq = None
        best_segment_validation_results = []
        candidate_read_found = False
        
        # collect candidate reads that overlap ref region of interest, then filter by size
        breakpoints = [pos for start, end in query.ref_segments for pos in (start, end)]
        bp_start, bp_stop = min(breakpoints), max(breakpoints)
        reads = self.bam_reader.candidate_read_seqs(
            query.chrom, bp_start, bp_stop, max_reads=self.max_reads_per_site)
        reads = [r for r in reads if len(r) >= len(query.sequence)]
        
        if reads:
            read_res = run_read_edlib(query.sequence, reads, self.read_error_threshold, self.min_read_support)
            candidate_read_found = read_res.was_read_found()  # TODO: also fold read_res.error into lowest_error
            if read_res.was_read_found() and read_res.error <= self.read_error_threshold:
                segment_validation_results = validate_segments_from_cigar(
                    Cigar.from_edlib(read_res.cigar), query.recon_segments,
                    error_threshold=self.read_error_threshold)
                if all(s.passed for s in segment_validation_results):
                    passed = True
                    if read_res.error < lowest_pass_error:
                        lowest_pass_error = read_res.error
                        best_segment_validation_results = segment_validation_results
                        validating_seq = read_res.matched_read_sequence
                if not passed:
                    best_segment_validation_results = segment_validation_results
        
        return QueryValidationInput(source=ValidationSource.READS, passed=passed,
                                    lowest_pass_error=lowest_pass_error, lowest_error=lowest_error,
                                    validating_seq=validating_seq, segment_results=best_segment_validation_results,
                                    overlapping_spanning_reads_found=bool(reads),
                                    candidate_read_found=candidate_read_found)


class AssemblyScorer(Scorer):
    """Validates a query against a FASTA (sample assembly, or the reference for the
    ambiguity check) via per-chromosome mappy alignments."""

    def __init__(self, config, chroms, fasta_path: str, forward_match_only: bool = False):
        super().__init__(config, chroms)
        self.forward_match_only = forward_match_only
        self.aligners = build_chrom_aligners(fasta_path, chroms, config.cache_dir)
        logger.info(f'Initialized per-chromosome aligners for {fasta_path}')

    def score_query(self, query: Query) -> QueryValidationInput:
        lowest_pass_error = 1.0
        lowest_error = 1.0
        passed = False
        validating_seq = None
        best_segment_validation_results = []
        best_strand_match = None

        aligner = self.aligners.get(query.chrom)
        if aligner is None:
            return QueryValidationInput(source=ValidationSource.ASSEMBLY)

        # gather and filter candidate alignments
        alignments: List[mappy.Alignment] = list(aligner.map(query.sequence)) # .map returns possible alignments, esp repetitive alignments?
        alignments = [a for a in alignments if abs(a.r_st - query.ref_start) <= self.location_tolerance]
        if self.forward_match_only:
            alignments = [a for a in alignments if a.strand == 1]

        # pick best candidate alignment based on overall and segment normalized match error
        for alignment in alignments:
            match_err = check_match(alignment, query.sequence)
            lowest_error = min(lowest_error, match_err)
            if match_err > self.match_error_threshold:
                continue

            segment_validation_results: List[SeqSegmentValidationResult] = \
                validate_segments_from_cigar(Cigar.from_mappy(alignment, len(query.sequence)),
                                             query.recon_segments,
                                             error_threshold=self.match_error_threshold)
            if all(s.passed for s in segment_validation_results):
                passed = True
                if match_err < lowest_pass_error:
                    lowest_pass_error = match_err
                    best_strand_match = alignment.strand   # record the best passing match's strand
                    best_segment_validation_results = segment_validation_results
                    matched = aligner.seq(alignment.ctg, alignment.r_st, alignment.r_en)
                    validating_seq = reverse_complement(matched) if alignment.strand == -1 else matched
            if not passed and match_err <= lowest_error:
                best_segment_validation_results = segment_validation_results

        return QueryValidationInput(
            source=ValidationSource.ASSEMBLY,
            passed=passed,
            lowest_pass_error=lowest_pass_error,
            lowest_error=lowest_error,
            validating_seq=validating_seq,
            segment_results=best_segment_validation_results,
            best_strand_match=best_strand_match,
        )


class EdlibScorer(Scorer):
    """Fallback for short queries mappy missed: expanding-window edlib search of a FASTA's bytes."""

    def __init__(self, config, chroms, fasta_path: str):
        super().__init__(config, chroms)
        self.min_edlib_query = config.min_edlib_query
        self.edlib_fallback_max_tolerance = config.edlib_fallback_max_tolerance
        self.sample_bytes = load_fasta_to_bytes(fasta_path, chroms)

    def score_query(self, query: Query) -> QueryValidationInput:
        if len(query.sequence) >= self.min_edlib_query:  # cost bound: the search is O(len * window)
            return QueryValidationInput(source=ValidationSource.EDLIB)

        lowest_pass_error = 1.0
        passed = False
        validating_seq = None
        segment_validation_results = []

        edlib_result = run_edlib_fallback(query.sequence, query.chrom, query.ref_start, query.buffer,
                                          self.edlib_fallback_max_tolerance, self.match_error_threshold,
                                          self.sample_bytes)
        if edlib_result is None:
            return QueryValidationInput(source=ValidationSource.EDLIB)

        if edlib_result.error <= self.match_error_threshold:
            segment_validation_results = validate_segments_from_cigar(
                Cigar.from_edlib(edlib_result.cigar), query.recon_segments,
                error_threshold=self.match_error_threshold)
            if all(s.passed for s in segment_validation_results):
                passed = True
                lowest_pass_error = edlib_result.error
                validating_seq = edlib_result.matched_target_sequence

        return QueryValidationInput(
            source=ValidationSource.EDLIB,
            passed=passed,
            lowest_pass_error=lowest_pass_error,
            lowest_error=edlib_result.error,
            validating_seq=validating_seq,
            segment_results=segment_validation_results,
        )
