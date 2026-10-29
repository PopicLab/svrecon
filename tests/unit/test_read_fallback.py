"""Tests ReadScorer when no read matches the allele within the error budget: the allele is measured
against a read's sequence between its flank anchors, or the query is inconclusive."""
import random
from types import SimpleNamespace

from svrecon.config import DEFAULTS
from svrecon.constants import QueryValidationReason, QueryValidationStatus
from svrecon.reconstruct import Query
from svrecon.scorers import reads
from svrecon.utils import reverse_complement

rng = random.Random(0)
random_seq = lambda n: ''.join(rng.choice('ACGT') for _ in range(n))

LEFT, A, B, C, RIGHT = (random_seq(n) for n in (100, 800, 600, 700, 100))
CLAIMED = LEFT + reverse_complement(A + B + C) + RIGHT  # the call: a plain inversion of ABC
# The sample: a rearrangement no INV matches, which also deletes the inner 38 bp of the right
# flank, as HG002 does at chr10:59,257,945.
TRUTH = LEFT + B + reverse_complement(A) + C + RIGHT[38:]
QUERY = Query(chrom='chrT', svtype='INV', grammar='~A~->~a~', sequence=CLAIMED, ref_start=1000,
              ref_end=1000 + len(CLAIMED), initial_segments=[],
              recon_segments=[(0, 100), (100, 2200), (2200, 2300)],
              ref_segments=[(1000, 1100), (1100, 3200), (3200, 3300)],
              ref_sequence=LEFT + A + B + C + RIGHT, buffer=100)


def score(monkeypatch, read_seqs):
    stub_reader = SimpleNamespace(candidate_read_seqs=lambda *args, **kwargs: read_seqs)
    monkeypatch.setattr(reads, 'BamReader', lambda bam_path: stub_reader)
    return reads.ReadScorer(SimpleNamespace(**{**DEFAULTS, 'bam': 'stub.bam'}), {'chrT'}).score_query(QUERY)


def test_miss_measured_between_flanks(monkeypatch):
    read = random_seq(3000) + reverse_complement(TRUTH) + random_seq(3000)  # sequenced from the other strand
    result = score(monkeypatch, [random_seq(9000), read])
    assert (result.status, result.reason) == (QueryValidationStatus.FAIL, QueryValidationReason.CIGAR_FAILED)
    assert result.best_matched_seq == TRUTH
    assert DEFAULTS['read_error_threshold'] < result.lowest_error < 1.0


def test_no_flanks_is_inconclusive(monkeypatch):
    result = score(monkeypatch, [random_seq(9000)])
    assert (result.status, result.reason) == (QueryValidationStatus.INCONCLUSIVE, QueryValidationReason.OTHER)
    assert result.best_matched_seq is None
