"""ReadScorer against the simulated BAM: one supporting read per SV.

TODO: test Cigar and the CigarValidations directly on hand-built CIGARs -- edlib picks
arbitrarily among equally optimal paths, so exact numbers can only be pinned without it.
"""
import pytest

from helpers import EXPECTED, SVID_IDS, SVIDS, get_reconstructed_query, get_reference_query
from svrecon.constants import QueryValidationReason, QueryValidationStatus, ValidationSource


@pytest.mark.parametrize('svid', SVIDS, ids=SVID_IDS)
def test_read_supports_reconstructed(read_scorer, reference, svs_by_svid, svid):
    """In: the alt allele. Out: PASS, zero error, one full-length match."""
    sv = svs_by_svid[svid]
    assert sv.svtype == EXPECTED[svid].svtype

    query = get_reconstructed_query(sv, reference)
    assert len(query) == EXPECTED[svid].reconstructed

    result = read_scorer.score_query(query)

    assert result.source is ValidationSource.READS
    assert result.status is QueryValidationStatus.PASS
    assert result.reason is QueryValidationReason.PASS
    assert result.passed is True
    assert result.lowest_error == 0.0
    assert result.lowest_pass_error == 0.0
    assert repr(result.cigar) == f'{len(query)}='


@pytest.mark.parametrize('svid', SVIDS, ids=SVID_IDS)
def test_read_rejects_reference(read_scorer, reference, svs_by_svid, svid):
    """In: the unrearranged reference. Out: FAIL -- the read carries the SV, not the reference."""
    query = get_reference_query(svs_by_svid[svid], reference)
    assert len(query) == EXPECTED[svid].reference

    result = read_scorer.score_query(query)

    assert result.source is ValidationSource.READS
    assert result.status is QueryValidationStatus.FAIL
    assert result.reason is QueryValidationReason.OTHER  # past edlib's bound: nothing aligned
    assert result.aligned is False
    assert result.passed is False
    assert result.lowest_error == 1.0
