"""Scorer contract: the Scorer base class and the QueryValidation result it returns."""
from dataclasses import dataclass, field
from typing import List, Optional

from svrecon.config import Config
from svrecon.constants import SubseqReason, SubseqStatus, ValidationSource
from svrecon.reconstruct import Query
from svrecon.scorers.utils import Cigar, SegmentValidation

# How much a status settles a query, least to most. A scorer that never ran reports the
# FAIL default, so it must rank below any scorer that actually looked at the sample.
_STATUS_RANK = {
    SubseqStatus.SKIPPED: 0,
    SubseqStatus.FAIL: 1,
    SubseqStatus.INCONCLUSIVE: 2,
    SubseqStatus.MATCH: 3,
    SubseqStatus.PASS: 4,
}


@dataclass
class QueryValidation:
    """One scorer's self-classified result."""
    source: ValidationSource
    passed: bool = False
    status: SubseqStatus = SubseqStatus.FAIL
    reason: SubseqReason = SubseqReason.OTHER
    lowest_pass_error: float = 1.0
    lowest_error: float = 1.0
    best_matched_seq: Optional[str] = None  # target sequence of the adopted alignment, kept on MATCH too
    segment_results: List[SegmentValidation] = field(default_factory=list)
    best_strand_match: Optional[int] = None  # alignmet only
    cigar: Optional[Cigar] = None  # the decisive alignment's CIGAR, for debugging

    def __gt__(self, other: 'QueryValidation') -> bool:
        """Ranks by status, then by lower error"""
        return ((_STATUS_RANK[self.status], -self.lowest_error)
                > (_STATUS_RANK[other.status], -other.lowest_error))

    def jsonify(self) -> dict:
        return {
            'source': self.source,
            'passed': self.passed,
            'status': self.status,
            'reason': self.reason,
            'cigar': self.cigar,
            'lowest_pass_error': self.lowest_pass_error,
            'lowest_error': self.lowest_error,
            'strand': self.best_strand_match,
            'segments': [s.jsonify() for s in self.segment_results],
        }


class Scorer:
    """Base for per-query validators; subclasses implement score_query(query) -> QueryValidation."""

    def __init__(self, config: Config, chroms: set[str]):
        # general shared params
        self.match_error_threshold = config.match_error_threshold # for both reference and assembly alignment
        self.location_tolerance = config.location_tolerance
        self.chroms = chroms
        # Sequence similarity always runs; the other two are opt-in via config.
        self.cigar_query_validations: List[CigarQueryValidation] = [
            SequenceSimilarityValidation(match_error_threshold=config.match_error_threshold),
        ]
        if config.max_indel_size is not None:
            self.cigar_query_validations.append(
                NoLargeIndelsValidation(max_size=config.max_indel_size))
        if config.junction_radius is not None:
            self.cigar_query_validations.append(
                JunctionValidation(radius=config.junction_radius,
                                   match_error_threshold=config.match_error_threshold))

    def score_query(self, query: Query) -> QueryValidation:
        raise NotImplementedError("Implement in subclass")

# TODO: look 
# Cigar Query Validations
@dataclass
class CigarQueryValidation:
    def validate(self, cigar: Cigar, query: Query) -> bool:
        raise NotImplementedError("Implement in Subclass")

@dataclass
class SequenceSimilarityValidation(CigarQueryValidation):
    """Bulk divergence over the whole query."""
    match_error_threshold: float

    def validate(self, cigar: Cigar, query: Query) -> bool:
        errors = cigar.total_error_bases + cigar.total_deleted_bases
        return errors / len(query) <= self.match_error_threshold

@dataclass
class NoLargeIndelsValidation(CigarQueryValidation):
    """Checks that no indel is size greater than ```max_size```"""
    max_size: int

    def validate(self, cigar: Cigar, query: Query) -> bool:
        indel_ops = {Cigar.INS} | Cigar.TARGET_ONLY_OPS
        return not any(length >= self.max_size
                       for op, length in cigar.cigartuples if op in indel_ops)

@dataclass
class JunctionValidation(CigarQueryValidation):
    """Local error around each breakpoint, where a wrong reconstruction shows up even when the
    allele as a whole aligns well."""
    radius: int
    match_error_threshold: float

    def validate(self, cigar: Cigar, query: Query) -> bool:
        breakpoints = {pos for start, end in query.recon_segments for pos in (start, end)}
        assert all(0 <= pos <= cigar.query_length for pos in breakpoints), \
            f'breakpoints {sorted(breakpoints)} outside the query [0,{cigar.query_length}]'

        for breakpoint_pos in sorted(breakpoints):
            if cigar.get_junction_error_rate(breakpoint_pos, self.radius) > self.match_error_threshold:
                return False
        return True