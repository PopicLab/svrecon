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
    CIGAR_FAILED = 'cigar_failed'
    OTHER = 'other'


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


class Tier(str, Enum):
    """Coarse evaluation tier for two-tier reporting."""
    READ = 'read'
    ASSEMBLY = 'assembly'
    MISS = 'miss'
