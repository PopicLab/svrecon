#!/usr/bin/env python3
"""Behavioral target for read-based SEGMENT validation (the junction -> segment refactor).

Black-box, no internals: build a complex SV (delINVdel), reconstruct its correct resulting
allele with `simulate_subsequences`, then push two "reads" through the real read path
(`run_read_edlib` -> CIGAR) and assert only the overall verdict of
`validate_segments_from_cigar`:

  * the CORRECT resulting subsequence                                   -> PASS
  * a WRONG subsequence that omits the SV's tiny right-flank deletion   -> FAIL

Why the wrong read matters: it differs from the correct allele by only the ~60 bp of the right
deletion, so its OVERALL edlib error is small (~3%) and a bulk/whole-sequence error check would
ACCEPT it. Rejecting it is the entire reason segment checking exists -- one segment (the right
del) is reconstructed wrong even though the sequence as a whole matches. `run_read_edlib` is
called with a large error_threshold so its bulk gate never pre-filters the wrong read;
`validate_segments_from_cigar` (at the real 0.1 threshold) is the sole arbiter.

`query.segments` is treated as opaque -- populated by `simulate_subsequences`, interpreted by
`validate_segments_from_cigar`. The test asserts pass/fail only, never segment structure.

This stays skipped until the refactor lands `svrecon.align.validate_segments_from_cigar`; once it
exists, the tests run (and will error, informatively, if `simulate_subsequences` isn't yet
populating `query.segments`).
"""
import random
import unittest

from svrecon.align import edlib_to_cigartuples, validate_segments_from_cigar
from svrecon.reads import run_read_edlib
from svrecon.reconstruct import simulate_subsequences
from svrecon.util import reverse_complement

# Seeded-random (NON-repetitive) reference: edlib is a real aligner, so a low-complexity/periodic
# reference lets it find degenerate alignments that absorb the wrong read's missing segment (making
# the wrong read look error-free). High complexity keeps alignments unambiguous. Long enough that
# the ~30 bp wrong segment is a small fraction of the allele -> overall error stays low while the
# segment check must still fail.
REF = ''.join(random.Random(20240730).choices('ACGT', k=2000))
BUF = 100_000  # >> len(REF): the single returned query is the whole reconstructed allele

# delINVdel over contiguous A / B / C in the middle of REF. Sizes chosen so that omitting C gives
# a small OVERALL error (~60/1840 = 3.2%, well under any bulk threshold) yet a large LOCAL error at
# the C junction (~60/300 = 20% in a 300 bp window) -- an unambiguous requirement for the segment
# check to fail the wrong read where a whole-sequence check would pass it.
A0, A1 = 800, 900      # left deletion (100 bp)
B0, B1 = 900, 1100     # inverted core (200 bp)
C0, C1 = 1100, 1160    # right deletion (60 bp) -- the "tiny right flank del"

READ_EDLIB_THRESHOLD = 1.0    # disable run_read_edlib's bulk gate; the segment check decides
SEGMENT_ERROR_THRESHOLD = 0.1


def _ref():
    return {'chr1': bytearray(REF, 'ascii')}


class _Rec:
    """Minimal pysam.VariantRecord stand-in; source region is 0-based [s, e)."""
    def __init__(self, s, e, op, svid='dS', svtype='delINVdel'):
        self.start = s + 1      # get_start_stop subtracts 1 -> s
        self.stop = e           # get_start_stop returns stop as-is -> e
        self.chrom = 'chr1'
        self.id = svid
        self.info = {'SVTYPE': svtype, 'SVID': svid, 'OP_TYPE': op}


def _reconstruct_query():
    """The correct resulting allele for the delINVdel (single full-allele query at BUF >> ref)."""
    recs = [_Rec(A0, A1, 'CUT'), _Rec(B0, B1, 'INV'), _Rec(C0, C1, 'CUT')]
    qs = simulate_subsequences(recs, BUF, _ref())
    seqs = {q.sequence for q in qs}
    assert len(seqs) == 1, f'expected one allele, got {len(seqs)}'
    return qs[0]


def _verdict(result):
    """Overall pass/fail from validate_segments_from_cigar. Mirrors the junction convention (a
    list of per-segment results each carrying `.passed`, all must pass); also tolerates a bare
    bool if the refactor returns one."""
    if isinstance(result, bool):
        return result
    return all(getattr(r, 'passed', r) for r in result)


class TestSegmentValidationReadPath(unittest.TestCase):
    def _read_passes(self, query, read_seq):
        res = run_read_edlib(query.sequence, [read_seq], error_threshold=READ_EDLIB_THRESHOLD)
        self.assertIsNotNone(res.cigar, 'read failed to align to the allele at all')
        cigartuples = edlib_to_cigartuples(res.cigar)
        results = validate_segments_from_cigar(cigartuples, query.segments,
                                               error_threshold=SEGMENT_ERROR_THRESHOLD)
        return _verdict(results)

    def test_correct_resulting_subsequence_passes(self):
        query = _reconstruct_query()
        # The exact reconstructed allele is, by definition, a correct read.
        self.assertTrue(self._read_passes(query, query.sequence))

    def test_wrong_subsequence_missing_tiny_right_del_fails(self):
        query = _reconstruct_query()
        # Correct allele deletes C; this read keeps C (right del not applied). It differs from the
        # allele by only ~60 bp -> small OVERALL error, but the right-del segment is wrong.
        wrong_read = REF[:A0] + reverse_complement(REF[B0:B1]) + REF[C0:]
        self.assertFalse(self._read_passes(query, wrong_read))

    def test_reference_read_without_sv_fails(self):
        query = _reconstruct_query()
        # A read carrying the reference (no SV at all): grossly wrong, must fail under any impl.
        self.assertFalse(self._read_passes(query, REF))


if __name__ == '__main__':
    unittest.main(verbosity=2)
