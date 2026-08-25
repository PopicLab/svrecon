"""Shared alignment primitives: normalized CIGAR representation, per-segment validation,
and edlib HW scoring."""
import re
from bisect import bisect_left, bisect_right
from dataclasses import dataclass
from typing import List, Tuple, Union

import edlib
import mappy
import pysam


class Cigar:
    """A CIGAR in normal form: pysam-ordered ``(op, length)`` tuples, forward-query
    orientation, 0 indexed, accounting for the full query (clipped flanks as explicit soft clips).
    Intervals are half-open [start, end)

    Range queries run on three prefix arrays sized by the op count, not the query length:
    ``query_starts[i]`` is the query index where op i begins (``query_starts[n]`` is the
    query length), ``error_prefix[i]`` / ``deletion_prefix[i]`` are the error / deleted
    bases in ops [0, i). Error ops are homogeneous, so a window that splits an op counts
    exactly its overlap.

    Official SAM/BAM CIGAR specification: https://samtools.github.io/hts-specs/SAMv1.pdf
    """

    # CIGAR operation codes (SAM/BAM spec), in code order.
    MATCH, INS, DEL, REF_SKIP, SOFT_CLIP, HARD_CLIP, PAD, SEQ_MATCH, SEQ_MISMATCH = range(9)
    QUERY_CONSUMING_OPS = {MATCH, INS, SOFT_CLIP, SEQ_MATCH, SEQ_MISMATCH}  # advance the query cursor
    TARGET_ONLY_OPS = {DEL, REF_SKIP}                                       # query gaps (deletions)
    ERROR_OPS = {INS, SOFT_CLIP, SEQ_MISMATCH}                              # query bases that aren't clean matches
    _EDLIB_OP_CODES = {'M': MATCH, '=': SEQ_MATCH, 'X': SEQ_MISMATCH, 'I': INS, 'D': DEL}
    _OP_CHARS = 'MIDNSHP=X'  # indexed by op code, for __repr__

    def __repr__(self) -> str:
        return ''.join(f'{length}{self._OP_CHARS[op]}' for op, length in self.cigartuples)

    def __init__(self, cigartuples: List[Tuple[int, int]]):
        self.cigartuples = cigartuples
        self.query_starts = [0]
        self.error_prefix = [0]
        self.deletion_prefix = [0]
        for op, length in cigartuples:
            self.query_starts.append(self.query_starts[-1] + (length if op in self.QUERY_CONSUMING_OPS else 0))
            self.error_prefix.append(self.error_prefix[-1] + (length if op in self.ERROR_OPS else 0))
            self.deletion_prefix.append(self.deletion_prefix[-1] + (length if op in self.TARGET_ONLY_OPS else 0))
        self.query_length = self.query_starts[-1]  # total query bases accounted for, clips included

    @property
    def total_error_bases(self) -> int:
        """Query bases that are clipped, inserted, or mismatched."""
        return self.error_prefix[-1]

    @property
    def total_deleted_bases(self) -> int:
        """Target bases the query skips over (deletions)."""
        return self.deletion_prefix[-1]

    def _errors_anchored(self, lo: int, hi: int) -> int:
        """Error bases among query bases [lo, hi), the end exclusive -- an error op occupies
        query positions, so only the bases inside the range count."""
        bounds = []
        for pos in (lo, hi):
            op_idx = bisect_right(self.query_starts, pos) - 1  # last op starting at or before pos
            if op_idx >= len(self.cigartuples):  # pos == query_length
                bounds.append(self.error_prefix[-1])
                continue
            # A split op contributes its overlap when it is an error op; ties on query_starts are
            # deletion ops, which consume no query and contribute 0 here.
            op = self.cigartuples[op_idx][0]
            partial = pos - self.query_starts[op_idx] if op in self.ERROR_OPS else 0
            bounds.append(self.error_prefix[op_idx] + partial)
        return bounds[1] - bounds[0]

    def _deletions_anchored(self, lo: int, hi: int) -> int:
        """Deleted target bases whose anchor -- the query base following the gap,
        ``query_starts[j]`` for deletion op j -- falls in [lo, hi], both ends inclusive: a
        deletion is a point between bases, so a gap sitting on either bound belongs to the
        range. ``hi == query_length`` therefore reaches a trailing deletion."""
        op_starts = self.query_starts[:-1]  # query_starts[j] is op j's anchor; drop the sentinel
        return self.deletion_prefix[bisect_right(op_starts, hi)] - self.deletion_prefix[bisect_left(op_starts, lo)]

    def get_junction_error_rate(self, point: int, radius: int) -> float:
        """Error rate over the window [point - radius, point + radius] clamped to the query:
        (error bases + deletions anchored in the window) / (window width + those deletions).
        Deletions at either edge count -- a junction's bounding gaps belong to it."""
        start = max(0, point - radius)
        end = min(self.query_length, point + radius)
        deletions = self._deletions_anchored(start, end)
        errors = self._errors_anchored(start, end) + deletions
        return errors / ((end - start) + deletions)

    @classmethod
    def from_mappy(cls, alignment: mappy.Alignment, query_len: int) -> 'Cigar':
        """Normalizes a mappy hit: flips (length, op) order, un-mirrors reverse-strand
        hits to forward-query order, adds the clips implicit in q_st/q_en.

        e.g. cigar=[(2000, M)], q_st=200, q_en=2200, query_len=2350
             -> [(S, 200), (M, 2000), (S, 150)]"""
        tuples = [(op, length) for length, op in alignment.cigar]
        if alignment.strand == -1:
            tuples.reverse()
        if alignment.q_st > 0:
            tuples.insert(0, (cls.SOFT_CLIP, alignment.q_st))
        if query_len - alignment.q_en > 0:
            tuples.append((cls.SOFT_CLIP, query_len - alignment.q_en))
        return cls(tuples)

    @classmethod
    def from_edlib(cls, cigar_str: str) -> 'Cigar':
        """Parses an edlib task='path' CIGAR string; HW mode is already full-query
        and forward, so parsing is the only normalization needed.

        e.g. '50=1X149=' -> [(=, 50), (X, 1), (=, 149)]"""
        return cls([(cls._EDLIB_OP_CODES[m.group(2)], int(m.group(1)))
                    for m in re.finditer(r'(\d+)([MIDX=])', cigar_str)])

    @classmethod
    def from_pysam(cls, read: pysam.AlignedSegment) -> 'Cigar':
        """Normalizes a pysam record: hard clips become soft clips (both mean unaligned
        original-read bases here); reverse-strand records flip to forward-read order.

        e.g. cigartuples=[(H, 100), (M, 50)], is_reverse=True -> [(M, 50), (S, 100)]"""
        tuples = [(cls.SOFT_CLIP if op == cls.HARD_CLIP else op, length) for op, length in read.cigartuples]
        if read.is_reverse:
            tuples.reverse()
        return cls(tuples)


@dataclass
class SegmentValidation:
    error: float
    passed: bool

    def jsonify(self) -> dict:
        return {'error': self.error, 'passed': self.passed}


def validate_segments_from_cigar(cigar: Cigar, segments: List[Tuple[int, int]],
                                 error_threshold: float = 0.1) -> List[SegmentValidation]:
    """Scores each segment [start, end) against the same error threshold; one result per segment."""
    rates = [cigar.get_window_error_rate(start, end) for start, end in segments]
    return [SegmentValidation(error=float(f'{rate:.4g}'), passed=rate <= error_threshold)
            for rate in rates]


@dataclass
class EdlibScoreResult:
    """One HW alignment of a query into a target (genomic window or read): normalized
    error rate, edlib CIGAR string, and the matched target substring."""
    error: float
    cigar: str  # edlib task='path' CIGAR string; parse with Cigar.from_edlib
    matched_target_sequence: str


def edlib_score(query_seq: str, target_seq: str, k: int = -1) -> Union[EdlibScoreResult, None]:
    """HW-align ``query_seq`` against a single ``target_seq`` and normalize to an
    error rate. Shared primitive for both assembly-window and read-based scoring;
    returns an ``EdlibScoreResult`` or ``None`` if no alignment was produced.

    ``k`` is edlib's max edit distance: alignments worse than ``k`` abort early
    and return ``None`` (edlib editDistance = -1). ``k=-1`` (default) is
    unbounded, preserving the assembly path's behavior; the read path passes a
    threshold-derived ``k`` so non-matching reads don't cost a full O(len*len)
    alignment -- the dominant cost when an SV is a read-mode miss."""
    result = edlib.align(query_seq.upper(), target_seq.upper(), mode="HW", task="path", k=k)
    if not result or result['editDistance'] < 0:
        return None

    match_start, match_end = result['locations'][0]  # best placement in target; edlib interval is inclusive
    if match_start is None:  # degenerate (empty query); task='path' otherwise always fills both
        return None

    matched_target_sequence = target_seq[match_start:match_end + 1]
    target_match_len = match_end - match_start + 1
    # TODO: revisit scoring. currently, normalize by the longer of query and matched span, so deletions widen thedenominator instead of inflating the rate.
    denominator = max(len(query_seq), target_match_len)
    error_rate = result['editDistance'] / denominator
    return EdlibScoreResult(error=error_rate, cigar=result['cigar'], matched_target_sequence=matched_target_sequence)
