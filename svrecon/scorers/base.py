"""Scorer contract: the Scorer base class and the QueryValidation result it returns."""
from dataclasses import dataclass, field
from typing import List, Optional

from svrecon.config import Config
from svrecon.constants import SubseqReason, SubseqStatus, ValidationSource
from svrecon.reconstruct import Query
from svrecon.scorers.utils import Cigar

# Cigar Query Validations
@dataclass
class CigarValidationResult:
    """One CigarQueryValidation's verdict. ``detail`` is free text so each check reports what it
    measured -- the error rate, the largest indel, every junction -- pass or fail alike."""
    passed: bool
    detail: str

    def jsonify(self) -> dict:
        return {'passed': self.passed, 'detail': self.detail}

@dataclass
class CigarQueryValidation:
    """Base for per-CIGAR checks; validate reports a verdict plus its own measurement."""

    def validate(self, cigar: Cigar, query: Query) -> CigarValidationResult:
        raise NotImplementedError("Implement in Subclass")

@dataclass
class CigarSequenceSimilarityValidation(CigarQueryValidation):
    """Bulk divergence over the whole query."""
    match_error_threshold: float

    def validate(self, cigar: Cigar, query: Query) -> CigarValidationResult:
        error = (cigar.total_error_bases + cigar.total_deleted_bases) / len(query)
        passed = error <= self.match_error_threshold
        return CigarValidationResult(
            passed, f'overall error {error:.4g} {"<=" if passed else ">"} {self.match_error_threshold}')

@dataclass
class CigarNoLargeIndelsValidation(CigarQueryValidation):
    """Checks that no indel is size greater than ```max_size```"""
    max_size: int

    def validate(self, cigar: Cigar, query: Query) -> CigarValidationResult:
        indel_ops = {Cigar.INS} | Cigar.TARGET_ONLY_OPS
        # max over (length, op) selects longest indel
        length, op = max(((length, op) for op, length in cigar.cigartuples if op in indel_ops), default=(0, None))
        passed = length < self.max_size
        kind = '' if op is None else (' insertion' if op == Cigar.INS else ' deletion')
        return CigarValidationResult(
            passed, f'largest indel {length} bp{kind} {"<" if passed else ">="} {self.max_size}')

@dataclass
class CigarJunctionValidation(CigarQueryValidation):
    """Local error around each breakpoint, where a wrong reconstruction shows up even when the
    allele as a whole aligns well."""
    radius: int
    match_error_threshold: float

    def validate(self, cigar: Cigar, query: Query) -> CigarValidationResult:
        breakpoints = {pos for start, end in query.recon_segments for pos in (start, end)}
        assert all(0 <= pos <= cigar.query_length for pos in breakpoints), \
            f'breakpoints {sorted(breakpoints)} outside the query [0,{cigar.query_length}]'
        if not breakpoints:
            return CigarValidationResult(True, 'no breakpoints')

        errors = [(pos, cigar.get_junction_error_rate(pos, self.radius)) for pos in sorted(breakpoints)]
        worst_pos, worst = max(errors, key=lambda pos_error: pos_error[-1])
        passed = worst <= self.match_error_threshold
        listing = ', '.join(f'{pos}:{error:.4g}' for pos, error in errors)

        return CigarValidationResult(
            passed, f'junction errors {{{listing}}}, worst {worst_pos}:{worst:.4g} '
                    f'{"<=" if passed else ">"} {self.match_error_threshold}')

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
    cigar_results: List[CigarValidationResult] = field(default_factory=list)
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
            'checks': [r.jsonify() for r in self.cigar_results],
        }


# How much a status settles a query, least to most. A scorer that never ran reports the
# FAIL default, so it must rank below any scorer that actually looked at the sample.
_STATUS_RANK = {
    SubseqStatus.SKIPPED: 0,
    SubseqStatus.FAIL: 1,
    SubseqStatus.INCONCLUSIVE: 2,
    SubseqStatus.MATCH: 3,
    SubseqStatus.PASS: 4,
}

class Scorer:
    """Base for per-query validators; subclasses implement score_query(query) -> QueryValidation."""

    def __init__(self, config: Config, chroms: set[str], error_threshold: Optional[float] = None):
        # general shared params
        self.match_error_threshold = config.match_error_threshold # for both reference and assembly alignment
        self.location_tolerance = config.location_tolerance
        self.chroms = chroms
        # Reads carry their own budget; every other scorer uses match_error_threshold.
        self.error_threshold = config.match_error_threshold if error_threshold is None else error_threshold
        # Sequence similarity always runs; the other two are opt-in via config.
        self.cigar_query_validations: List[CigarQueryValidation] = [
            CigarSequenceSimilarityValidation(match_error_threshold=self.error_threshold),
        ]
        if config.max_indel_size is not None:
            self.cigar_query_validations.append(
                CigarNoLargeIndelsValidation(max_size=config.max_indel_size))
        if config.junction_radius is not None:
            self.cigar_query_validations.append(
                CigarJunctionValidation(radius=config.junction_radius,
                                   match_error_threshold=self.error_threshold))

    def score_cigar(self, cigar: Cigar, query: Query) -> List[CigarValidationResult]:
        """Every configured check's verdict for one alignment, in configured order."""
        return [validation.validate(cigar, query) for validation in self.cigar_query_validations]

    def score_query(self, query: Query) -> QueryValidation:
        raise NotImplementedError("Implement in subclass")

