"""Shared alignment primitives: normalized CIGAR representation, per-segment validation,
and edlib HW scoring."""
import re
from bisect import bisect_right
from collections import deque
from dataclasses import dataclass
from typing import List, Tuple, Union

import edlib
import mappy
import pysam


class Cigar:
    """A CIGAR in normal form: pysam-ordered ``(op, length)`` tuples, forward-query orientation,
    0 indexed, full query length (clipped flanks as explicit soft clips).

    Official SAM/BAM CIGAR specification: https://samtools.github.io/hts-specs/SAMv1.pdf
    """

    # CIGAR operation codes (SAM/BAM spec), in code order.
    MATCH, INS, DEL, REF_SKIP, SOFT_CLIP, HARD_CLIP, PAD, SEQ_MATCH, SEQ_MISMATCH = range(9)
    QUERY_CONSUMING_OPS = {MATCH, INS, SOFT_CLIP, SEQ_MATCH, SEQ_MISMATCH}  # advances the query cursor
    TARGET_ONLY_OPS = {DEL, REF_SKIP}                                       # query gaps
    ERROR_OPS = {INS, SOFT_CLIP, SEQ_MISMATCH}                              # query bases that aren't clean matches
    COLUMN_CONSUMING_OPS = QUERY_CONSUMING_OPS | TARGET_ONLY_OPS            # advances the column cursor
    COLUMN_ERROR_OPS = ERROR_OPS | TARGET_ONLY_OPS                          # error columns: query-side errors + deletions
    _EDLIB_OP_CODES = {'M': MATCH, '=': SEQ_MATCH, 'X': SEQ_MISMATCH, 'I': INS, 'D': DEL}
    _OP_CHARS = 'MIDNSHP=X'  # cigar ops, indexed by op code

    def __repr__(self) -> str:
        return ''.join(f'{length}{self._OP_CHARS[op]}' for op, length in self.cigartuples)

    def __init__(self, cigartuples: List[Tuple[int, int]]):
        self.cigartuples = cigartuples
        # attributes following mappy.Alignment, of aligned block
        self.blen = sum(op_len for op, op_len in cigartuples  # aligned block: M/I/D/=/X, clips excluded
                        if op in self.COLUMN_CONSUMING_OPS and op != self.SOFT_CLIP)
        self.NM = sum(op_len for op, op_len in cigartuples    # edit distance: mismatches + ins + del
                      if op in {self.SEQ_MISMATCH, self.INS, self.DEL})
        # attributes from aligning entire sequences
        self.num_columns_total = sum(op_len for op, op_len in self.cigartuples if op in self.COLUMN_CONSUMING_OPS)
        self.num_columns_error = sum(op_len for op, op_len in self.cigartuples if op in self.COLUMN_ERROR_OPS)

    @property
    def block_aligned_error_rate(self) -> float:
        """Divergence within the aligned block: NM / blen. Clips don't count."""
        return self.NM / self.blen

    @property
    def whole_query_error_rate(self) -> float:
        """Full alignment divergence: every error column over every column, clips included."""
        return self.num_columns_error / self.num_columns_total

    def get_left_clip(self) -> int:
        """Length of the leading soft clip"""
        op, length = self.cigartuples[0]
        return length if op == self.SOFT_CLIP else 0

    def get_right_clip(self) -> int:
        """Length of the trailing soft clip"""
        op, length = self.cigartuples[-1]
        return length if op == self.SOFT_CLIP else 0

    def get_error_windows(self, window_size: int, error: int) -> List[Tuple[int, int]]:
        """return all window ofs ```window_size``` alignment columns across the CIGAR"""
        windows = []
        cigar_operations = [[op, op_capacity] for op, op_capacity in self.cigartuples if op in self.COLUMN_CONSUMING_OPS]
        window: deque = deque()                 
        num_window_columns = num_window_errors = 0
        query_pos = 0                           
        cigar_idx = 0                           

        # build the initial window, consuming all we can of each operation
        while num_window_columns < window_size and cigar_idx < len(cigar_operations):
            op, op_capacity = cigar_operations[cigar_idx]
            take = min(op_capacity, window_size - num_window_columns)
            window.append([op, take])
            num_window_columns += take
            num_window_errors += take if op in self.COLUMN_ERROR_OPS else 0
            cigar_operations[cigar_idx][1] -= take
            if cigar_operations[cigar_idx][1] == 0:
                cigar_idx += 1
        windows.append((query_pos, num_window_errors))
        if num_window_columns < window_size:        # fewer columns than one window
            return windows

        # traverse cigar: drop the window's trailing piece, refill from the remaining capacity from the right
        while cigar_idx < len(cigar_operations):
            op, num_dropped_ops = window.popleft()
            num_window_columns -= num_dropped_ops
            num_window_errors -= num_dropped_ops if op in self.COLUMN_ERROR_OPS else 0
            query_pos += num_dropped_ops if op in self.QUERY_CONSUMING_OPS else 0

            while num_window_columns < window_size and cigar_idx < len(cigar_operations):
                op, op_capacity = cigar_operations[cigar_idx]
                take = min(op_capacity, window_size - num_window_columns)
                window.append([op, take])
                num_window_columns += take
                num_window_errors += take if op in self.COLUMN_ERROR_OPS else 0
                cigar_operations[cigar_idx][1] -= take
                if cigar_operations[cigar_idx][1] == 0:
                    cigar_idx += 1
            if num_window_columns < window_size:    # out of columns; no further full window exists
                return windows
            windows.append((query_pos, num_window_errors))
        return windows

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
        """Parses an edlib CIGAR string.

        e.g. '50=1X149=' -> [(=, 50), (X, 1), (=, 149)]"""
        return cls([(cls._EDLIB_OP_CODES[m.group(2)], int(m.group(1)))
                    for m in re.finditer(r'(\d+)([MIDX=])', cigar_str)])


@dataclass
class EdlibScoreResult:
    """One HW alignment of a query into a target (genomic window or read): normalized
    error rate, the parsed CIGAR, and the matched target substring."""
    error: float
    cigar: Cigar  # parsed from edlib's task='path' CIGAR; error is its NM/blen
    matched_target_sequence: str


def edlib_score(query_seq: str, target_seq: str, mode="HW", k: int = -1) -> Union[EdlibScoreResult, None]:
    """HW-align ``query_seq`` against a single ``target_seq`` and normalize to an
    error rate. Shared primitive for both assembly-window and read-based scoring;
    returns an ``EdlibScoreResult`` or ``None`` if no alignment was produced.

    ``k`` is edlib's max edit distance: alignments worse than ``k`` abort early
    and return ``None`` (edlib editDistance = -1). ``k=-1`` (default) is
    unbounded, preserving the assembly path's behavior; the read path passes a
    threshold-derived ``k`` so non-matching reads don't cost a full O(len*len)
    alignment -- the dominant cost when an SV is a read-mode miss."""
    result = edlib.align(query_seq.upper(), target_seq.upper(), mode=mode, task="path", k=k)
    if not result or result['editDistance'] < 0:
        return None

    match_start, match_end = result['locations'][0]  # best placement in target; edlib interval is inclusive
    if match_start is None:  # degenerate (empty query); task='path' otherwise always fills both
        return None

    matched_target_sequence = target_seq[match_start:match_end + 1]

    cigar = Cigar.from_edlib(result['cigar'])
    return EdlibScoreResult(error=cigar.whole_query_error_rate, cigar=cigar,
                            matched_target_sequence=matched_target_sequence)
