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

    def score_query(self, query: Query) -> QueryValidation:
        raise NotImplementedError("Implement in subclass")
