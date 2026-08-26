"""Scorer contract: the Scorer base class and the QueryValidation result it returns."""
from dataclasses import dataclass, field
from typing import List, Optional, Tuple

from svrecon.config import Config
from svrecon.constants import QueryValidationReason, QueryValidationStatus, ValidationSource
from svrecon.reconstruct import Query
from svrecon.scorers.utils import Cigar

PASS = QueryValidationStatus.PASS
FAIL = QueryValidationStatus.FAIL
INCONCLUSIVE = QueryValidationStatus.INCONCLUSIVE

MAPPABLE_BASES = 'ACGT'  # every other IUPAC code (N above all) cannot be aligned against

# Cigar Query Validations
@dataclass
class CigarQueryValidationResult:
    """One CigarQueryValidation's verdict. ``detail`` is free text so each check reports what it
    measured -- the error rate, the largest indel, every junction -- whatever the status."""
    status: QueryValidationStatus
    detail: str

    def jsonify(self) -> dict:
        return {'status': self.status, 'detail': self.detail}

@dataclass
class CigarQueryValidation:
    """Base for per-CIGAR checks; validate reports a status plus its own measurement."""

    def validate(self, cigar: Cigar, query: Query) -> CigarQueryValidationResult:
        raise NotImplementedError("Implement in Subclass")

@dataclass
class CigarQuerySequenceSimilarityValidation(CigarQueryValidation):
    """Bulk divergence over the whole query."""
    match_error_threshold: float

    def validate(self, cigar: Cigar, query: Query) -> CigarQueryValidationResult:
        error = (cigar.total_error_bases + cigar.total_deleted_bases) / len(query)
        passed = error <= self.match_error_threshold
        return CigarQueryValidationResult(
            PASS if passed else FAIL,
            f'overall error {error:.4g} {"<=" if passed else ">"} {self.match_error_threshold}')

@dataclass
class CigarQueryNoLargeIndelsValidation(CigarQueryValidation):
    """Checks that no indel is size greater than ```max_size```"""
    max_size: int

    def validate(self, cigar: Cigar, query: Query) -> CigarQueryValidationResult:
        indel_ops = {Cigar.INS} | Cigar.TARGET_ONLY_OPS
        # max over (length, op) selects longest indel
        length, op = max(((length, op) for op, length in cigar.cigartuples if op in indel_ops), default=(0, None))
        passed = length < self.max_size
        kind = '' if op is None else (' insertion' if op == Cigar.INS else ' deletion')
        return CigarQueryValidationResult(
            PASS if passed else FAIL,
            f'largest indel {length} bp{kind} {"<" if passed else ">="} {self.max_size}')

@dataclass
class CigarQueryJunctionValidation(CigarQueryValidation):
    """Local error around each breakpoint, where a wrong reconstruction shows up even when the
    allele as a whole aligns well."""
    radius: int
    match_error_threshold: float

    def validate(self, cigar: Cigar, query: Query) -> CigarQueryValidationResult:
        breakpoints = {pos for start, end in query.recon_segments for pos in (start, end)}
        assert all(0 <= pos <= cigar.query_length for pos in breakpoints), \
            f'breakpoints {sorted(breakpoints)} outside the query [0,{cigar.query_length}]'
        breakpoints -= {0, cigar.query_length} # exclude start and end of query in the check
        if not breakpoints:
            return CigarQueryValidationResult(PASS, 'no breakpoints')

        errors = [(pos, cigar.get_junction_error_rate(pos, self.radius)) for pos in sorted(breakpoints)]
        worst_pos, worst = max(errors, key=lambda pos_error: pos_error[-1])
        passed = worst <= self.match_error_threshold
        listing = ', '.join(f'{pos}:{error:.4g}' for pos, error in errors)

        return CigarQueryValidationResult(
            PASS if passed else FAIL,
            f'junction errors {{{listing}}}, worst {worst_pos}:{worst:.4g} '
            f'{"<=" if passed else ">"} {self.match_error_threshold}')

@dataclass
class CigarQueryMappableValidation(CigarQueryValidation):
    """Inconclusive if fraction of mappable bases of the query (A/C/G/T) is not above ```min_mappable_fraction```"""
    min_mappable_fraction: float

    def validate(self, cigar: Cigar, query: Query) -> CigarQueryValidationResult:
        sequence = query.sequence.upper()
        # str.count runs in C; a per-base Python loop would rerun for every candidate alignment
        mappable = sum(sequence.count(base) for base in MAPPABLE_BASES)
        fraction = mappable / len(sequence) if sequence else 0.0
        passed = fraction >= self.min_mappable_fraction
        return CigarQueryValidationResult(
            PASS if passed else INCONCLUSIVE,
            f'mappable {fraction:.4g} {">=" if passed else "<"} {self.min_mappable_fraction}')

@dataclass
class QueryValidation:
    """One scorer's self-classified result."""
    source: ValidationSource
    passed: bool = False
    status: QueryValidationStatus = QueryValidationStatus.FAIL
    reason: QueryValidationReason = QueryValidationReason.OTHER
    lowest_pass_error: float = 1.0
    lowest_error: float = 1.0
    best_matched_seq: Optional[str] = None  # target sequence of the adopted alignment, kept on an aligned fail too
    cigar_results: List[CigarQueryValidationResult] = field(default_factory=list)
    best_strand_match: Optional[int] = None  # alignmet only
    cigar: Optional[Cigar] = None  # the decisive alignment's CIGAR, for debugging

    @property
    def aligned(self) -> bool:
        """Whether an alignment was found -- a CIGAR failure is one a check then rejected."""
        return self.passed or self.reason is QueryValidationReason.CIGAR_FAILED

    @property
    def rank(self) -> int:
        """How much this result settles the query, least to most. The FAIL/OTHER default a scorer
        that never ran reports is lowest; an aligned fail contradicts, so it beats inconclusive."""
        if self.status is QueryValidationStatus.PASS:
            return 3
        if self.status is QueryValidationStatus.FAIL and self.aligned:
            return 2
        if self.status is QueryValidationStatus.INCONCLUSIVE:
            return 1
        return 0

    def __gt__(self, other: 'QueryValidation') -> bool:
        """Ranks by how much the result settles the query, then by lower error"""
        return (self.rank, -self.lowest_error) > (other.rank, -other.lowest_error)

    def jsonify(self) -> dict:
        return {
            'source': self.source,
            'passed': self.passed,
            'status': self.status,
            'reason': self.reason,
            'aligned': self.aligned,
            'cigar': self.cigar,
            'lowest_pass_error': self.lowest_pass_error,
            'lowest_error': self.lowest_error,
            'strand': self.best_strand_match,
            'checks': [r.jsonify() for r in self.cigar_results],
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
        # Sequence similarity always runs; the rest are opt-in via config.
        self.cigar_query_validations: List[CigarQueryValidation] = [
            CigarQuerySequenceSimilarityValidation(match_error_threshold=self.error_threshold),
        ]
        if config.max_indel_size is not None:
            self.cigar_query_validations.append(
                CigarQueryNoLargeIndelsValidation(max_size=config.max_indel_size))
        if config.junction_radius is not None:
            self.cigar_query_validations.append(
                CigarQueryJunctionValidation(radius=config.junction_radius,
                                   match_error_threshold=self.error_threshold))
        if config.min_mappable_fraction is not None:
            self.cigar_query_validations.append(
                CigarQueryMappableValidation(min_mappable_fraction=config.min_mappable_fraction))

    def score_cigar(self, cigar: Cigar,
                    query: Query) -> Tuple[QueryValidationStatus, List[CigarQueryValidationResult]]:
        """score cigar and query against all validation methods. PASS if all checks pass. INCONCLUSIVE
        if at least one check is INCONCLUSIVE. else Fail """
        results = [validation.validate(cigar, query) for validation in self.cigar_query_validations]
        statuses = {result.status for result in results}
        if INCONCLUSIVE in statuses:
            return INCONCLUSIVE, results
        return (FAIL if FAIL in statuses else PASS), results

    def score_query(self, query: Query) -> QueryValidation:
        raise NotImplementedError("Implement in subclass")

