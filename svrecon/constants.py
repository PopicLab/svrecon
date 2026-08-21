"""Shared enums for SV scoring classification (status/reason/outcome)."""
from enum import Enum


class SubseqStatus(str, Enum):
    """Per-subsequence verdict."""
    PASS = 'pass'
    MATCH = 'match'  # aligned somewhere within the error threshold, but a segment failed
    FAIL = 'fail'
    INCONCLUSIVE = 'inconclusive'
    SKIPPED = 'skipped'  # neither --sample nor --bam -- reconstructed but deliberately not validated


class SubseqReason(str, Enum):
    """Why a subsequence received its SubseqStatus."""
    PASS = 'pass'
    REFERENCE_MATCH = 'reference_match'
    SEGMENT_FAILED = 'segment_failed'
    OVER_ERROR_THRESHOLD = 'over_error_threshold'
    INCONCLUSIVE = 'inconclusive'
    SKIPPED = 'skipped'  # neither --sample nor --bam
    OTHER = 'other'


class Outcome(str, Enum):
    """SV-level roll-up outcome across all of its reconstructed subsequences."""
    HIT = 'hit'
    MISS = 'miss'
    INCONCLUSIVE = 'inconclusive'
    SKIPPED = 'skipped'  # neither --sample nor --bam -- reconstructed but deliberately not validated


class ValidationSource(str, Enum):
    """Which tier ultimately validated (or was attempted for) a subsequence."""
    READS = 'reads'
    ASSEMBLY = 'assembly'
    EDLIB = 'edlib'


class Tier(str, Enum):
    """Coarse evaluation tier for two-tier reporting."""
    READ = 'read'
    ASSEMBLY = 'assembly'
    MISS = 'miss'
