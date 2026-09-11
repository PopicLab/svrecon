"""Scorer contract: the Scorer base class and the QueryValidation result it returns."""
import re
from dataclasses import dataclass, field
from typing import List, Optional, Tuple

from svrecon.config import Config
from svrecon.constants import QueryValidationReason, QueryValidationStatus, SOURCE_RANK, ValidationSource
from svrecon.reconstruct import Query
from svrecon.scorers.utils import Cigar

PASS = QueryValidationStatus.PASS
FAIL = QueryValidationStatus.FAIL
INCONCLUSIVE = QueryValidationStatus.INCONCLUSIVE

MAPPABLE_BASES = 'ACGT'  # every other IUPAC code (N above all) cannot be aligned against

# Cigar Query Checks
@dataclass
class CigarQueryCheckResult:
    """One CigarQueryCheck's verdict. ``detail`` is free text so each check reports what it
    measured -- the error rate, the largest indel, every junction -- whatever the status."""
    status: QueryValidationStatus
    detail: str

    def jsonify(self) -> dict:
        return {'status': self.status, 'detail': self.detail}

@dataclass
class CigarQueryCheck:
    """Base for per-CIGAR checks; validate reports a status plus its own measurement."""

    def validate(self, cigar: Cigar, query: Query) -> CigarQueryCheckResult:
        raise NotImplementedError("Implement in Subclass")

@dataclass
class CigarQueryAlignmentSimilarityCheck(CigarQueryCheck):
    """Bulk divergence over the whole query."""
    match_error_threshold: float

    def validate(self, cigar: Cigar, query: Query) -> CigarQueryCheckResult:
        error = cigar.NM / cigar.blen
        passed = error <= self.match_error_threshold
        return CigarQueryCheckResult(
            PASS if passed else FAIL,
            f'alignment error {error:.4g} {"<=" if passed else ">"} {self.match_error_threshold}')

@dataclass
class CigarQueryNoLargeErrorsCheck(CigarQueryCheck):
    """Checks that no error run -- indel or soft-clipped flank -- is ``max_size`` or longer."""
    max_size: int

    def validate(self, cigar: Cigar, query: Query) -> CigarQueryCheckResult:
        indel_ops = {Cigar.INS} | Cigar.TARGET_ONLY_OPS
        large_errors = [(length, op) for op, length in cigar.cigartuples
                        if op in indel_ops and length >= self.max_size]
        large_errors += [(clip, Cigar.SOFT_CLIP) for clip in (cigar.get_left_clip(), cigar.get_right_clip())
                         if clip >= self.max_size]
        if large_errors:
            large_error_details = ', '.join(f'{length}{Cigar._OP_CHARS[op]}' for length, op in large_errors)
            return CigarQueryCheckResult(FAIL, f'error runs {{{large_error_details}}} >= {self.max_size}')
        return CigarQueryCheckResult(PASS, f'no error run >= {self.max_size}')

@dataclass
class CigarQueryMaximumErrorWindowCheck(CigarQueryCheck):
    """Checks that no window of ``window_size`` alignment columns holds ``error`` or more."""
    max_window_size: int
    max_window_error: int

    def validate(self, cigar: Cigar, query: Query) -> CigarQueryCheckResult:
        windows = cigar.get_error_windows(self.max_window_size, self.max_window_error)
        offending = [(pos, errors) for pos, errors in windows if errors >= self.max_window_error]
        if offending:
            listing = ', '.join(f'{pos}:{errors}' for pos, errors in offending)
            return CigarQueryCheckResult(
                FAIL, f'error windows {{{listing}}} >= {self.max_window_error} per {self.max_window_size}')
        worst_pos, worst_errors = max(windows, key=lambda window: window[1], default=(0, 0))
        return CigarQueryCheckResult(
            PASS, f'no window >= {self.max_window_error} errors per {self.max_window_size} '
                  f'(worst {worst_pos}:{worst_errors})')

@dataclass
class CigarQueryMappableCheck(CigarQueryCheck):
    """
    pass if every configured condition is met 
    - fraction of mappable bases of the query (A/C/G/T) is above ```min_mappable_fraction```
    - no large contiguous nonmappable region of at least length ```max_unmappable_size```
    otherwise inconclusive
    """
    min_mappable_fraction: Optional[float]
    max_unmappable_size: Optional[int]

    def validate(self, cigar: Cigar, query: Query) -> CigarQueryCheckResult:
        sequence = query.sequence.upper()
        checks = []  # (passed, detail) pairs, one per condition that's actually configured

        if self.min_mappable_fraction is not None:
            mappable = sum(sequence.count(base) for base in MAPPABLE_BASES)
            mappable_fraction = mappable / len(sequence) if sequence else 0.0
            mapple_fraction_pass = mappable_fraction >= self.min_mappable_fraction
            checks.append((mapple_fraction_pass, f'mappable {mappable_fraction:.4g} '
                          f'{">=" if mapple_fraction_pass else "<"} {self.min_mappable_fraction}'))

        if self.max_unmappable_size is not None:
            longest_gap = max((len(m.group()) for m in re.finditer(f'[^{MAPPABLE_BASES}]+', sequence)), default=0)
            gap_pass = longest_gap < self.max_unmappable_size
            checks.append((gap_pass, f'longest nonmappable run {longest_gap} '
                          f'{"<" if gap_pass else ">="} {self.max_unmappable_size}'))

        passed = all(check_pass for check_pass, _ in checks)
        return CigarQueryCheckResult(PASS if passed else INCONCLUSIVE, ', '.join(d for _, d in checks))

@dataclass
class QueryValidation:
    """One scorer's self-classified result."""
    source: ValidationSource
    passed: bool = False
    status: QueryValidationStatus = QueryValidationStatus.FAIL
    reason: QueryValidationReason = QueryValidationReason.NOT_ATTEMPTED
    lowest_pass_error: float = 1.0
    lowest_error: float = 1.0
    best_matched_seq: Optional[str] = None  # target sequence of the adopted alignment, kept on an aligned fail too
    cigar_results: List[CigarQueryCheckResult] = field(default_factory=list)
    best_strand_match: Optional[int] = None  # alignmet only
    cigar: Optional[Cigar] = None  # the decisive alignment's CIGAR, for debugging

    @property
    def aligned(self) -> bool:
        """Whether an alignment was found -- a CIGAR failure is one a check then rejected."""
        return self.best_matched_seq is not None

    @property
    def status_rank(self) -> int:
        if self.status is QueryValidationStatus.PASS:
            return 3
        elif self.status is QueryValidationStatus.FAIL and self.aligned:
            return 2
        elif self.status is QueryValidationStatus.INCONCLUSIVE:
            return 1
        elif self.status is QueryValidationStatus.FAIL and not self.aligned:
            return 0
        else:
            raise RuntimeError(f"Unexpected branch, {self.status=} is not set")

    @property
    def source_rank(self) -> int:
        return SOURCE_RANK[self.source]

    def __gt__(self, other: 'QueryValidation') -> bool:
        """Ranks by how much the result settles the query, then by validation source tier,
        then by lower error"""
        return ((self.status_rank, self.source_rank, -self.lowest_error) >
               (other.status_rank, other.source_rank, -other.lowest_error))

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
        self.cigar_query_checks: List[CigarQueryCheck] = [
            CigarQueryAlignmentSimilarityCheck(match_error_threshold=self.error_threshold),
        ]
        if config.max_contiguous_error is not None:
            self.cigar_query_checks.append(
                CigarQueryNoLargeErrorsCheck(max_size=config.max_contiguous_error))
        if (config.max_window_size is None) != (config.max_window_error is None):
            raise ValueError('max_window_size and max_window_error must be set together '
                             f'(got {config.max_window_size!r} and {config.max_window_error!r})')
        if config.max_window_size is not None:
            self.cigar_query_checks.append(
                CigarQueryMaximumErrorWindowCheck(max_window_size=config.max_window_size,
                                                  max_window_error=config.max_window_error))
        if config.min_mappable_fraction is not None or config.max_unmappable_size is not None:
            self.cigar_query_checks.append(
                CigarQueryMappableCheck(min_mappable_fraction=config.min_mappable_fraction,
                                        max_unmappable_size=config.max_unmappable_size))

    def score_cigar(self, cigar: Cigar,
                    query: Query) -> Tuple[QueryValidationStatus, List[CigarQueryCheckResult]]:
        """score cigar and query against all checks. PASS if all checks pass. INCONCLUSIVE
        if at least one check is INCONCLUSIVE. else Fail """
        results = [check.validate(cigar, query) for check in self.cigar_query_checks]
        statuses = {result.status for result in results}
        if INCONCLUSIVE in statuses:
            return INCONCLUSIVE, results
        return (FAIL if FAIL in statuses else PASS), results

    def score_query(self, query: Query) -> QueryValidation:
        raise NotImplementedError("Implement in subclass")
