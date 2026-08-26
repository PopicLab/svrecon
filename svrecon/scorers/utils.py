"""Shared alignment primitives: normalized CIGAR representation, per-segment validation,
and edlib HW scoring."""
import re
from bisect import bisect_right
from dataclasses import dataclass
from typing import List, Tuple, Union

import edlib
import mappy
import pysam


class Cigar:
    """A CIGAR in normal form: pysam-ordered ``(op, length)`` tuples, forward-query orientation,
    0 indexed, full query length (clipped flanks as explicit soft clips). Intervals are half-open
    [start, end).

    Two coordinate systems, one per prefix-array pair below, sized by op count (not length):

    - query position: index into the query sequence. Only query-consuming ops advance it, so a
      deletion (target-only) has no query position of its own.
    - alignment column: index into the alignment itself. Query-consuming and target-only ops
      both advance it, so a deletion occupies columns despite consuming no query.

    Each array below is indexed by op index i in ``cigartuples``. Worked example throughout:
    ``2S3=2D1X4=`` (ops 0-4: soft clip, match, deletion, mismatch, match).

    - ``query_starts[i]`` -- the query position where op i begins (``query_starts[n]`` is the
      total query length).
        ex. ``query_starts = [0,2,5,5,6,10]``
    - ``column_starts[i]`` -- the column where op i begins (``column_starts[n]`` is the total
      column count).
        ex. ``column_starts = [0,2,5,7,8,12]``
    - ``error_prefix[i]`` -- error bases accumulated over ops [0, i).
        ex. ``error_prefix = [0,2,2,2,3,3]``
    - ``deletion_prefix[i]`` -- deleted (target-only) bases accumulated over ops [0, i).
        ex. ``deletion_prefix = [0,0,0,2,2,2]``

    Ops are homogeneous, so a window that splits one counts exactly its overlap.

    Official SAM/BAM CIGAR specification: https://samtools.github.io/hts-specs/SAMv1.pdf
    """

    # CIGAR operation codes (SAM/BAM spec), in code order.
    MATCH, INS, DEL, REF_SKIP, SOFT_CLIP, HARD_CLIP, PAD, SEQ_MATCH, SEQ_MISMATCH = range(9)
    QUERY_CONSUMING_OPS = {MATCH, INS, SOFT_CLIP, SEQ_MATCH, SEQ_MISMATCH}  # advance the query cursor
    TARGET_ONLY_OPS = {DEL, REF_SKIP}                                       # query gaps (deletions)
    ERROR_OPS = {INS, SOFT_CLIP, SEQ_MISMATCH}                              # query bases that aren't clean matches
    COLUMN_CONSUMING_OPS = QUERY_CONSUMING_OPS | TARGET_ONLY_OPS            # advance the column cursor
    _EDLIB_OP_CODES = {'M': MATCH, '=': SEQ_MATCH, 'X': SEQ_MISMATCH, 'I': INS, 'D': DEL}
    _OP_CHARS = 'MIDNSHP=X'  # indexed by op code, for __repr__

    def __repr__(self) -> str:
        return ''.join(f'{length}{self._OP_CHARS[op]}' for op, length in self.cigartuples)

    def __init__(self, cigartuples: List[Tuple[int, int]]):
        self.cigartuples = cigartuples
        self.query_starts = [0]
        self.column_starts = [0]
        self.error_prefix = [0]
        self.deletion_prefix = [0]
        for op, length in cigartuples:
            self.query_starts.append(self.query_starts[-1] + (length if op in self.QUERY_CONSUMING_OPS else 0))
            self.column_starts.append(self.column_starts[-1] + (length if op in self.COLUMN_CONSUMING_OPS else 0))
            self.error_prefix.append(self.error_prefix[-1] + (length if op in self.ERROR_OPS else 0))
            self.deletion_prefix.append(self.deletion_prefix[-1] + (length if op in self.TARGET_ONLY_OPS else 0))
        self.query_length = self.query_starts[-1]  # total query bases accounted for, clips included
        self.column_length = self.column_starts[-1]  # total alignment columns, query bases and gaps alike

    @property
    def total_error_bases(self) -> int:
        """Query bases that are clipped, inserted, or mismatched."""
        return self.error_prefix[-1]

    @property
    def total_deleted_bases(self) -> int:
        """Target bases the query skips over (deletions)."""
        return self.deletion_prefix[-1]

    @property
    def error_rate(self) -> float:
        """Bulk error over the whole query: (error bases + deletions) / query length."""
        return (self.total_error_bases + self.total_deleted_bases) / self.query_length

    def _column_of(self, query_pos: int) -> int:
        """The alignment column occupied by query base ``query_pos``. bisect_right (not _left)
        skips past any zero-query-width deletion ops tied at the same query_starts value, landing
        on the op that actually holds the base."""
        op_idx = bisect_right(self.query_starts, query_pos) - 1
        return self.column_starts[op_idx] + (query_pos - self.query_starts[op_idx])

    def _at_column(self, column: int) -> Tuple[int, int]:
        """``(error columns, query columns)`` in the half-open range [0, column). A column is an
        error if it's a deletion or an error base -- anything that isn't a clean match."""
        op_idx = bisect_right(self.column_starts, column) - 1
        if op_idx >= len(self.cigartuples):  # column == column_length
            return self.error_prefix[-1] + self.deletion_prefix[-1], self.query_starts[-1]
        op = self.cigartuples[op_idx][0]
        partial = column - self.column_starts[op_idx]  # a split op contributes its overlap
        errors = self.error_prefix[op_idx] + self.deletion_prefix[op_idx]
        if op in self.ERROR_OPS or op in self.TARGET_ONLY_OPS:
            errors += partial
        queries = self.query_starts[op_idx] + (partial if op in self.QUERY_CONSUMING_OPS else 0)
        return errors, queries

    def get_junction_error_rate(self, point: int, radius: int) -> float:
        """Error rate around query base ``point``: from its alignment column, expand ``radius``
        columns each way (closed, clamped to the alignment), score (error columns) / (query bases)
        in that span.

        e.g. ``3D5=`` (columns ``DDD=====``), point=0 (the first ``=``), radius=1 spans ``D==``:
        one deletion column reached out of three, over two query bases -- 1/2, not 3/2 or 0."""
        assert 0 <= point < self.query_length, f'point {point} outside the query [0,{self.query_length})'
        center = self._column_of(point)
        start = max(0, center - radius)
        end = min(self.column_length - 1, center + radius)  # closed: the column at `end` counts
        errors_end, queries_end = self._at_column(end + 1)
        errors_start, queries_start = self._at_column(start)
        # `center` is a query-consuming column inside the window, so the denominator is never 0
        return (errors_end - errors_start) / (queries_end - queries_start)

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
    # TODO: revisit scoring. currently, normalize by the longer of query and matched span, so deletions widen thedenominator instead of inflating the rate.
    denominator = len(query_seq)
    error_rate = result['editDistance'] / denominator
    return EdlibScoreResult(error=error_rate, cigar=result['cigar'], matched_target_sequence=matched_target_sequence)
