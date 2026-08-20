"""Shared scorer utilities: normalized CIGAR representation."""
import re
from functools import cached_property
from typing import List, Tuple, TYPE_CHECKING

import mappy
import numpy as np
import pysam



class Cigar:
    """A CIGAR in normal form: pysam-ordered ``(op, length)`` tuples, forward-query
    orientation, 0 indexed, accounting for the full query (clipped flanks as explicit soft clips).
    Intervals are half-open [start, end)

    Official SAM/BAM CIGAR specification: https://samtools.github.io/hts-specs/SAMv1.pdf
    """

    # CIGAR operation codes (SAM/BAM spec), in code order.
    MATCH, INS, DEL, REF_SKIP, SOFT_CLIP, HARD_CLIP, PAD, SEQ_MATCH, SEQ_MISMATCH = range(9)
    QUERY_CONSUMING_OPS = {MATCH, INS, SOFT_CLIP, SEQ_MATCH, SEQ_MISMATCH}  # advance the query cursor
    TARGET_ONLY_OPS = {DEL, REF_SKIP}                                       # query gaps (deletions)
    ERROR_OPS = {INS, SOFT_CLIP, SEQ_MISMATCH}                              # query bases that aren't clean matches
    _EDLIB_OP_CODES = {'M': MATCH, '=': SEQ_MATCH, 'X': SEQ_MISMATCH, 'I': INS, 'D': DEL}

    def __init__(self, cigartuples: List[Tuple[int, int]]):
        self.cigartuples = cigartuples

    @cached_property
    def query_length(self) -> int:
        """Total query bases the CIGAR accounts for, clips included."""
        return sum(length for op, length in self.cigartuples if op in self.QUERY_CONSUMING_OPS)

    @cached_property
    def error_mask(self) -> np.ndarray:
        """boolean bit mask. mask[i] = 1 iff query base i is clipped, inserted, or mismatched."""
        mask = np.zeros(self.query_length, dtype=np.uint8)
        pos = 0
        for op, length in self.cigartuples:
            if op in self.QUERY_CONSUMING_OPS:
                if op in self.ERROR_OPS:
                    mask[pos:pos + length] = 1
                pos += length
        return mask

    @cached_property
    def deletion_lengths(self) -> np.ndarray:
        """int32, length query_length + 1. dels[i] = n indicates a deletion of n target bases
        preceding query base i."""
        dels = np.zeros(self.query_length + 1, dtype=np.int32)
        pos = 0
        for op, length in self.cigartuples:
            if op in self.TARGET_ONLY_OPS:
                dels[pos] += length
            elif op in self.QUERY_CONSUMING_OPS:
                pos += length
        return dels

    def get_window_error_rate(self, start: int, end: int) -> float:
        """Error rate over [start, end):
        (error bases + deletions anchored in [start, end)) / ((end - start) + those deletions)."""
        del_len = int(self.deletion_lengths[start:end].sum())
        errors = int(self.error_mask[start:end].sum()) + del_len
        return errors / ((end - start) + del_len) # TODO: why do we include this?

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
