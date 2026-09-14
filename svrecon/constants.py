"""Shared enums for SV scoring classification (status/reason/outcome)."""
from enum import Enum


class QueryValidationStatus(str, Enum):
    """Per-query verdict."""
    PASS = 'pass'
    FAIL = 'fail'
    INCONCLUSIVE = 'inconclusive'


class QueryValidationReason(str, Enum):
    """Why a query received its QueryValidationStatus."""
    PASS = 'pass'
    CIGAR_FAILED = 'cigar failed'
    CIGAR_INCONCLUSIVE = 'cigar inconclusive'
    NOT_ATTEMPTED = 'not attempted'  # no scorer ran: not configured, or declined (e.g. a cost bound)
    REFERENCE_AMBIGUOUS = 'reference ambiguous'  # the allele also matches the reference

    # Read specific
    NO_SPANNING_READS = 'no spanning reads'
    NO_PASSING_READ = 'all reads are bellow error threshold, so we score against the first read found'
    NO_MATCHED_BASE_PAIRS_IN_FAILING_READS = 'all reads are below the error threshold, so we score against the first read found. However, the first read has no matching bases within the reference span'

    # Assembly Specific
    NO_MATCHING_CONTIG = 'no matching contig'

    # Edlib specific
    NO_ALIGNMENT_FOUND = 'no alignment found'  # the expanding-window search found nothing within max_tolerance


class Outcome(str, Enum):
    """SV-level roll-up outcome across all of its reconstructed queries."""
    HIT = 'hit'
    MISS = 'miss'
    INCONCLUSIVE = 'inconclusive'  # untestable: no scorer reached a verdict, or none was configured


class ValidationSource(str, Enum):
    """Which tier ultimately validated (or was attempted for) a query."""
    READS = 'reads'
    ASSEMBLY = 'assembly'
    EDLIB = 'edlib'


SOURCE_RANK = {ValidationSource.EDLIB: 0, ValidationSource.ASSEMBLY: 1, ValidationSource.READS: 2}


class Tier(str, Enum):
    """Coarse evaluation tier for two-tier reporting."""
    READ = 'read'
    ASSEMBLY = 'assembly'
    MISS = 'miss'
