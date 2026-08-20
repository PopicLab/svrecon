"""Scorer contract: the Scorer base class and the QueryValidationInput result it returns."""
from dataclasses import dataclass, field
from typing import List, Optional

from svrecon.config import Config
from svrecon.constants import SubseqReason, SubseqStatus, ValidationSource
from svrecon.reconstruct import Query
from svrecon.scorers.utils import SegmentValidation


@dataclass
class QueryValidationInput:
    """One scorer's self-classified result."""
    source: ValidationSource
    passed: bool = False
    status: SubseqStatus = SubseqStatus.FAIL
    reason: SubseqReason = SubseqReason.OTHER
    lowest_pass_error: float = 1.0
    lowest_error: float = 1.0
    validating_seq: Optional[str] = None
    segment_results: List[SegmentValidation] = field(default_factory=list)
    best_strand_match: Optional[int] = None  # alignmet only

    def jsonify(self) -> dict:
        return {
            'source': self.source,
            'passed': self.passed,
            'status': self.status,
            'reason': self.reason,
            'lowest_pass_error': self.lowest_pass_error,
            'lowest_error': self.lowest_error,
            'strand': self.best_strand_match,
            'segments': [s.jsonify() for s in self.segment_results],
        }


class Scorer:
    """Base for per-query validators; subclasses implement score_query(query) -> QueryValidationInput."""

    def __init__(self, config: Config, chroms: set[str]):
        # general shared params
        self.match_error_threshold = config.match_error_threshold # for both reference and assembly alignment
        self.location_tolerance = config.location_tolerance
        self.chroms = chroms

    def score_query(self, query: Query) -> QueryValidationInput:
        raise NotImplementedError("Implement in subclass")
