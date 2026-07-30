"""Shared enums for SV scoring classification (status/reason/outcome)."""
from enum import Enum


class SubseqStatus(str, Enum):
    """Per-subsequence verdict, exactly one of these three."""
    PASS = 'pass'
    FAIL = 'fail'
    INCONCLUSIVE = 'inconclusive'


class SubseqReason(str, Enum):
    """Why a subsequence received its SubseqStatus."""
    PASS = 'pass'
    REFERENCE_MATCH = 'reference_match'
    JUNCTION_FAILED = 'junction_failed'
    OVER_ERROR_THRESHOLD = 'over_error_threshold'
    NO_READS = 'no_reads'
    INCONCLUSIVE = 'inconclusive'
    NO_ALIGNER = 'no_aligner'
    NO_ALIGNMENT_IN_WINDOW = 'no_alignment_in_window'


class Outcome(str, Enum):
    """SV-level roll-up outcome across all of its reconstructed subsequences."""
    HIT = 'hit'
    MISS = 'miss'
    INCONCLUSIVE = 'inconclusive'
