#!/usr/bin/env python3
"""Unit tests for the reconstruction scorer's pure helpers.

Run:  python -m pytest tests/      (from the repo root)
  or: python -m unittest discover tests

Scope: the functions that can be exercised without VCF/BAM/aligner fixtures --
`reverse_complement` (svrecon.util), `edlib_score`/`validate_segments_from_cigar`
(svrecon.align) and `run_read_edlib` (svrecon.reads). Several cases are regression tests for bugs fixed while
building read-based evaluation; those carry a "regression" note. The reason/tier
classification inside `score_sv` is intentionally not covered here -- it is
inlined and would require full VCF/BAM/aligner fixtures.
"""
import unittest

from svrecon.util import reverse_complement
from svrecon.align import edlib_score, validate_segments_from_cigar
from svrecon.reads import run_read_edlib, ReadEdlibResult

# BAM CIGAR op codes used to build fixtures below.
SEQ_MATCH, SEQ_MISMATCH = 7, 8


class TestReverseComplement(unittest.TestCase):
    """Single consolidated revcomp (fast str.translate); accepts str or a list of
    single-char strings and always returns a str."""

    def test_str_input(self):
        self.assertEqual(reverse_complement('AACG'), 'CGTT')

    def test_list_input_returns_str(self):
        # Callers pass a list slice of the reference (list of single chars).
        rc = reverse_complement(list('AACG'))
        self.assertIsInstance(rc, str)
        self.assertEqual(rc, 'CGTT')

    def test_lowercase_and_n(self):
        self.assertEqual(reverse_complement('acgtN'), 'Nacgt')

    def test_double_revcomp_is_identity(self):
        seq = 'ACGTTGCACCGGATN'
        self.assertEqual(reverse_complement(reverse_complement(seq)), seq)

    def test_slice_assignment_contract(self):
        # score_alignments.py inverts a region via `buf[a:b] = reverse_complement(...)`;
        # assigning a str to a list slice must splice into per-char elements.
        buf = list('AAAAAA')
        buf[1:5] = reverse_complement(list('AACG'))
        self.assertEqual(buf, list('ACGTTA'))

    def test_clip_concat_contract(self):
        # Inverted dispersion clips are wrapped in list() so they concatenate with
        # the list-based haplotype sequence.
        clip = list(reverse_complement(list('AACG')))
        self.assertEqual(['X'] + clip + ['Y'], ['X', 'C', 'G', 'T', 'T', 'Y'])


class TestEdlibScore(unittest.TestCase):
    """Shared HW (infix) edlib primitive returning an EdlibScoreResult or None."""

    def test_exact_infix_match(self):
        res = edlib_score('ACGTACGT', 'GGGACGTACGTGGG')
        self.assertIsNotNone(res)
        self.assertEqual(res.error, 0.0)
        self.assertIsNotNone(res.cigar)

    def test_single_mismatch_has_small_error(self):
        res = edlib_score('ACGTACGT', 'ACGTTCGT')  # one substitution
        self.assertIsNotNone(res)
        self.assertTrue(0.0 < res.error < 0.2)

    def test_k_bound_aborts_to_none(self):
        # editDistance beyond k -> edlib reports -1 -> None (the bound that keeps
        # read-mode misses from costing a full O(len*len) alignment).
        self.assertIsNone(edlib_score('ACGTACGTAC', 'TTTTTTTTTT', k=1))


class TestRunReadEdlib(unittest.TestCase):
    ALT = 'A' * 40 + 'C' * 20  # non-palindromic; revcomp is clearly different

    def test_returns_dataclass(self):
        out = run_read_edlib(self.ALT, [self.ALT], 0.1)
        self.assertIsInstance(out, ReadEdlibResult)

    def test_forward_exact_match(self):
        res = run_read_edlib(self.ALT, [self.ALT], 0.1)
        self.assertIsNotNone(res.cigar)
        self.assertEqual(res.error, 0.0)
        self.assertEqual(res.n_tried, 1)

    def test_reverse_strand_read_matches(self):
        # Regression (sv1427): a read sequenced from the opposite strand is the
        # reverse complement of the reference-oriented alt. Forward-only alignment
        # scored it ~1.0 and the call missed; run_read_edlib must try both strands.
        rc_read = reverse_complement(self.ALT)
        self.assertEqual(edlib_score(self.ALT, rc_read).error, 1.0)  # forward alone: miss
        res = run_read_edlib(self.ALT, [rc_read], 0.1)
        self.assertIsNotNone(res.cigar)
        self.assertLessEqual(res.error, 0.1)
        self.assertEqual(res.n_tried, 1)

    def test_short_read_is_length_filtered(self):
        # A read too short to host the whole alt can't contain it; it's skipped
        # before alignment, so n_tried stays 0 (distinct from "aligned but failed").
        res = run_read_edlib(self.ALT, ['A' * 20], 0.1)
        self.assertIsNone(res.cigar)
        self.assertEqual(res.n_tried, 0)

    def test_min_support_controls_early_exit(self):
        # min_support only gates the early return, not whether a result comes back.
        res1 = run_read_edlib(self.ALT, [self.ALT, self.ALT], 0.1, min_support=1)
        res2 = run_read_edlib(self.ALT, [self.ALT, self.ALT], 0.1, min_support=2)
        self.assertEqual(res1.n_tried, 1)
        self.assertEqual(res2.n_tried, 2)


class TestValidateSegmentsFromCigar(unittest.TestCase):
    """The validator returns an ordered list of per-segment SeqSegmentValidationResult objects
    (not a bool). Overall verdict is all(s.passed for s in results); an empty list (no in-scope
    segment) passes. It short-circuits on the first failing segment. A zero-width segment (i, i)
    is the junction-style point case, checked over the window [i - radius, i + radius]; a
    (start, end) segment extends the window to [start - radius, end + radius]."""

    def test_clean_match_passes(self):
        # 200 exact matches; a point segment mid-alignment sees zero local error.
        results = validate_segments_from_cigar([(200, SEQ_MATCH)], [(100, 100)], radius=50, error_threshold=0.1)
        self.assertEqual(len(results), 1)
        self.assertTrue(results[0].passed)
        self.assertEqual(results[0].error, 0.0)  # native number
        self.assertTrue(all(s.passed for s in results))

    def test_mismatch_cluster_at_segment_fails(self):
        # 20 mismatches inside a 100bp window (err 0.2 > 0.1) -> the point segment fails.
        cig = [(50, SEQ_MATCH), (20, SEQ_MISMATCH), (130, SEQ_MATCH)]
        results = validate_segments_from_cigar(cig, [(60, 60)], radius=50, error_threshold=0.1)
        self.assertEqual(len(results), 1)
        self.assertFalse(results[0].passed)
        self.assertAlmostEqual(results[0].error, 0.2, places=6)  # native number
        self.assertFalse(all(s.passed for s in results))

    def test_range_segment_pools_error_over_window(self):
        # A real (non-point) segment [50,150) with radius 25 -> window [25,175], 150 bases; the
        # 30 mismatches inside it give 30/150 = 0.2 > 0.1 -> the segment fails. Exercises the
        # range window that a point junction could not.
        cig = [(60, SEQ_MATCH), (30, SEQ_MISMATCH), (110, SEQ_MATCH)]
        results = validate_segments_from_cigar(cig, [(50, 150)], radius=25, error_threshold=0.1)
        self.assertEqual(len(results), 1)
        self.assertFalse(results[0].passed)
        self.assertAlmostEqual(results[0].error, 30 / 150, places=6)

    def test_short_circuits_on_first_failure(self):
        # Two segments, the first fails -> the list ends at it; the second is not checked.
        cig = [(50, SEQ_MATCH), (20, SEQ_MISMATCH), (130, SEQ_MATCH)]
        results = validate_segments_from_cigar(cig, [(60, 60), (500, 500)], radius=50, error_threshold=0.1)
        self.assertEqual(len(results), 1)
        self.assertFalse(results[0].passed)

    def test_out_of_scope_segment_skipped(self):
        # A segment beyond this alignment span is not counted (empty -> passes).
        results = validate_segments_from_cigar([(200, SEQ_MATCH)], [(10000, 10000)], radius=50)
        self.assertEqual(results, [])
        self.assertTrue(all(s.passed for s in results))

    def test_no_segments_returns_empty_and_passes(self):
        results = validate_segments_from_cigar([(200, SEQ_MATCH)], [], radius=50)
        self.assertEqual(results, [])
        self.assertTrue(all(s.passed for s in results))


if __name__ == '__main__':
    unittest.main(verbosity=2)
