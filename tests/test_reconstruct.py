#!/usr/bin/env python3
"""Unit tests for `simulate_subsequences` -- the alt-allele reconstruction.

Run:  python -m unittest discover tests   (from the repo root)

Strategy: use a BUFFER larger than the reference so `get_changed_subsequences`
returns the whole reconstructed allele as a single query. Then the reconstructed
`sequence` can be compared directly against an independently hand-built expected
allele (string ops on the reference), and junction offsets are absolute positions
in that allele -- no window arithmetic. Each fundamental operation
(DEL/INV/COPY-PASTE/COPYinv-PASTE/CUT-PASTE/CUTinv-PASTE), a couple of multi-op
orderings, and the compound dupINVdup grammar are exercised.

Coordinate contract (verified empirically): a record's source region is 0-based
half-open ``[s, e)`` passed as ``start=s+1`` (1-based POS) / ``stop=e``; with a
buffer >= len(ref) the working ``offset`` clamps to 0, so op coordinates and
insert ``TARGET`` are absolute 0-based indices into the reference.
"""
import unittest

from svrecon.reconstruct import simulate_subsequences
from svrecon.util import reverse_complement

# Deterministic, non-homopolymer reference so slices and their revcomps are distinct.
_B = 'ACGT'
REF_S = ''.join(_B[(i * 7 + i // 4) % 4] for i in range(200))
BUF = 10_000  # >> len(REF_S): the single returned query is the full allele


def ref():
    return {'chr1': bytearray(REF_S, 'ascii')}


class _Rec:
    """Minimal pysam.VariantRecord stand-in. Source region is 0-based [s, e)."""
    def __init__(self, s, e, op, svtype, target=None, svid='sv1', insord=None):
        self.start = s + 1      # get_start_stop subtracts 1 -> s
        self.stop = e           # get_start_stop returns stop as-is -> e
        self.chrom = 'chr1'
        self.id = svid
        self.info = {'SVTYPE': svtype, 'SVID': svid, 'OP_TYPE': op}
        if target is not None:
            self.info['TARGET'] = target
        if insord is not None:
            self.info['INSORD'] = insord


def one(recs):
    """Run the reconstruction and return one query. With BUF >> len(ref) every
    changed interval's window spans the whole allele, so multiple changed regions
    (e.g. cut+paste) yield identical full-allele queries; assert they agree and
    return the first."""
    qs = simulate_subsequences(recs, BUF, ref())
    seqs = {q['sequence'] for q in qs}
    assert len(seqs) == 1, f'queries disagree on the allele ({len(seqs)} distinct): {qs}'
    return qs[0]


class TestFundamentalOps(unittest.TestCase):
    def test_deletion(self):
        q = one([_Rec(60, 80, 'CUT', 'DEL')])
        self.assertEqual(q['sequence'], REF_S[:60] + REF_S[80:])

    def test_inversion(self):
        q = one([_Rec(60, 80, 'INV', 'INV')])
        self.assertEqual(q['sequence'], REF_S[:60] + reverse_complement(REF_S[60:80]) + REF_S[80:])
        # INV marks the first and last base of the inverted block.
        self.assertEqual(q['junctions'], [60, 79])

    def test_tandem_dup(self):
        # DUP sets TARGET=rec.stop -> copy inserted right after the source (tandem).
        q = one([_Rec(60, 80, 'COPY-PASTE', 'DUP')])
        self.assertEqual(q['sequence'], REF_S[:80] + REF_S[60:80] + REF_S[80:])
        self.assertEqual(q['junctions'], [80, 99])  # endpoints of the inserted copy

    def test_dispersed_copy_paste(self):
        # Copy [60,80) forward to index 120 (downstream); source stays in place.
        q = one([_Rec(60, 80, 'COPY-PASTE', 'dDUP', target=120)])
        self.assertEqual(q['sequence'], REF_S[:120] + REF_S[60:80] + REF_S[120:])
        self.assertEqual(q['junctions'], [120, 139])

    def test_dispersed_copyinv_paste(self):
        # Copy reverse-complement of [60,80) to index 120.
        q = one([_Rec(60, 80, 'COPYinv-PASTE', 'INV_dDUP', target=120)])
        self.assertEqual(q['sequence'], REF_S[:120] + reverse_complement(REF_S[60:80]) + REF_S[120:])
        self.assertEqual(q['junctions'], [120, 139])

    def test_cut_paste(self):
        # Cut [60,80) and paste forward at index 120.
        q = one([_Rec(60, 80, 'CUT-PASTE', 'nrTRA', target=120)])
        # source removed; copy inserted before original index 120 (indices unaffected: 120>80).
        self.assertEqual(q['sequence'], REF_S[:60] + REF_S[80:120] + REF_S[60:80] + REF_S[120:])

    def test_cutinv_paste(self):
        q = one([_Rec(60, 80, 'CUTinv-PASTE', 'INV_nrTRA', target=120)])
        self.assertEqual(q['sequence'], REF_S[:60] + REF_S[80:120] + reverse_complement(REF_S[60:80]) + REF_S[120:])


class TestOperationOrdering(unittest.TestCase):
    def test_two_dispersed_inserts_are_order_independent(self):
        # Two forward copies to different targets; sorted-descending processing must
        # keep both insertion points correct regardless of index shifts.
        r1 = _Rec(20, 40, 'COPY-PASTE', 'dDUP', target=90, svid='svX', insord=0)
        r2 = _Rec(50, 60, 'COPY-PASTE', 'dDUP', target=150, svid='svX', insord=1)
        q = one([r1, r2])
        # Apply right-to-left so earlier splices don't shift later (lower) targets.
        expected = REF_S[:150] + REF_S[50:60] + REF_S[150:]
        expected = expected[:90] + REF_S[20:40] + expected[90:]
        self.assertEqual(q['sequence'], expected)

    def test_inversion_plus_downstream_dup(self):
        inv = _Rec(60, 80, 'INV', 'INV', svid='svY')
        dup = _Rec(100, 110, 'COPY-PASTE', 'dDUP', target=140, svid='svY')
        q = one([inv, dup])
        expected = REF_S[:60] + reverse_complement(REF_S[60:80]) + REF_S[80:]
        expected = expected[:140] + REF_S[100:110] + expected[140:]
        self.assertEqual(q['sequence'], expected)


class TestDupINVdupGrammar(unittest.TestCase):
    """dupINVdup grammar ABC -> A c b a C, where c=RC(C), b=RC(B), a=RC(A). The two
    NOVEL adjacencies (the only non-reference joins) are A|c and a|C -- exactly what
    the junction check must validate. The three fragments that build it: copy C
    (inverted) to the A|B boundary, invert B in place, copy A (inverted) to the B|C
    boundary."""

    def _fragments(self, cCa_target, aA_target, b_lo, b_hi):
        # cCa_target: where RC(C) is inserted; aA_target: where RC(A) is inserted;
        # [b_lo,b_hi): the in-place inverted B.
        return [
            _Rec(130, 170, 'COPYinv-PASTE', 'dupINVdup', target=cCa_target, svid='dS'),  # C -> c
            _Rec(b_lo, b_hi, 'INV', 'dupINVdup', svid='dS'),                              # B -> b
            _Rec(40, 70, 'COPYinv-PASTE', 'dupINVdup', target=aA_target, svid='dS'),      # A -> a
        ]

    def test_engine_reconstructs_clean_acbac(self):
        # Contiguous A=[40,70) B=[70,130) C=[130,170); c before 70, a before 130.
        q = one(self._fragments(cCa_target=70, aA_target=130, b_lo=70, b_hi=130))
        A, B, C = REF_S[40:70], REF_S[70:130], REF_S[130:170]
        expected = (REF_S[:40] + A + reverse_complement(C) + reverse_complement(B)
                    + reverse_complement(A) + C + REF_S[170:])
        self.assertEqual(q['sequence'], expected)
        # Novel junctions land at the A|c and a|C boundaries.
        self.assertIn(70, q['junctions'])
        self.assertIn(70 + len(C), q['junctions'])  # start of a-block (a|C is at its far end)

    @unittest.expectedFailure
    def test_groovi_fragment_coords_leave_stray_junction_bases(self):
        # Documents a groovi bug (seq/stitch_rules.py dupINVdup): B=[p1_1+1, p2_0-1)
        # drops the two breakpoint bases and insert targets are asymmetric
        # (pos2_0-1 vs pos1_1), so the allele is A c [stray] b a [stray] C, not A c b a C.
        # Reproduce groovi's EXACT coords and assert the CLEAN allele -> currently fails.
        p1_1, p2_0 = 70, 130
        q = one(self._fragments(cCa_target=p1_1, aA_target=p2_0 - 1, b_lo=p1_1 + 1, b_hi=p2_0 - 1))
        A, B, C = REF_S[40:70], REF_S[70:130], REF_S[130:170]
        expected = (REF_S[:40] + A + reverse_complement(C) + reverse_complement(B)
                    + reverse_complement(A) + C + REF_S[170:])
        self.assertEqual(q['sequence'], expected)


if __name__ == '__main__':
    unittest.main(verbosity=2)
