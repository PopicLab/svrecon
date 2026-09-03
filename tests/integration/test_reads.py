"""ReadScorer against the simulated BAM: one supporting read per SV."""
import pytest

from constants import CHROM, QUERY_BUFFER_SIZE, READS_BAM
from helpers import (EXPECTED_SVTYPES, SIMULATED_SVIDS, SVID_TEST_IDS, get_true_reconstructed_query,
                     get_reference_query, get_scorer_config)
from svrecon.constants import QueryValidationReason, QueryValidationStatus, ValidationSource
from svrecon.scorers.reads import ReadScorer


@pytest.fixture(scope='module')
def read_scorer() -> ReadScorer:
    """ReadScorer over BAM"""
    return ReadScorer(get_scorer_config(bam=str(READS_BAM)), chroms={CHROM})


@pytest.mark.parametrize('svid', SIMULATED_SVIDS, ids=SVID_TEST_IDS)
def test_read_supports_reconstructed(read_scorer, assembly_sequence,
                                     svs_by_svid, svid):
    """In: the alt allele. Out: PASS, zero error, one full-length match."""
    sv = svs_by_svid[svid]
    assert sv.svtype == EXPECTED_SVTYPES[svid]

    query = get_true_reconstructed_query(sv, assembly_sequence)

    result = read_scorer.score_query(query)

    assert result.source is ValidationSource.READS
    assert result.status is QueryValidationStatus.PASS
    assert result.reason is QueryValidationReason.PASS
    assert result.passed is True
    assert result.lowest_error == 0.0
    assert result.lowest_pass_error == 0.0


@pytest.mark.parametrize('svid', SIMULATED_SVIDS, ids=SVID_TEST_IDS)
def test_read_rejects_reference(read_scorer, reference_sequence,
                                svs_by_svid, svid):
    """In: the unrearranged reference. Out: FAIL -- the read carries the SV, not the reference."""
    query = get_reference_query(svs_by_svid[svid], reference_sequence)

    result = read_scorer.score_query(query)

    assert result.source is ValidationSource.READS
    assert result.status is QueryValidationStatus.FAIL
    assert result.reason is QueryValidationReason.OTHER  # past edlib's bound: nothing aligned
    assert result.aligned is False
    assert result.passed is False
