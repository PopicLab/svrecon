"""Scorer contract: the Scorer base class and the QueryValidationInput result it returns."""
from dataclasses import dataclass, field
from typing import List, Optional

from svrecon.config import Config
from svrecon.constants import SubseqReason, SubseqStatus, ValidationSource
from svrecon.reconstruct import Query
from svrecon.scorers.utils import SegmentValidation


@dataclass
class QueryValidationInput:
    """One scorer's self-classified result for one query; each scorer sets its own
    status/reason, so the folding in QueryValidation stays scorer-agnostic."""
    source: ValidationSource
    passed: bool = False
    status: SubseqStatus = SubseqStatus.FAIL
    reason: SubseqReason = SubseqReason.OTHER
    lowest_pass_error: float = 1.0
    lowest_error: float = 1.0
    validating_seq: Optional[str] = None
    segment_results: List[SegmentValidation] = field(default_factory=list)
    best_strand_match: Optional[int] = None  # assembly only


class Scorer:
    """Base for per-query validators; subclasses implement score_query(query) -> QueryValidationInput."""

    def __init__(self, config: Config, chroms: set[str]):
        # general shared params
        self.match_error_threshold = config.match_error_threshold # for both reference and assembly alignment
        self.location_tolerance = config.location_tolerance
        self.chroms = chroms

    def score_query(self, query: Query) -> QueryValidationInput:
        raise NotImplementedError("Implement in subclass")
