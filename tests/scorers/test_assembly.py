"""AssemblyScorer against the assembly insilicoSV built: mappy alignment of each allele."""
import pytest

from helpers import EXPECTED, SVID_IDS, SVIDS, get_reconstructed_query, get_reference_query
from svrecon.constants import QueryValidationReason, QueryValidationStatus, ValidationSource


@pytest.mark.parametrize('svid', SVIDS, ids=SVID_IDS)
def test_alignment_supports_reconstructed(assembly_scorer, reference, svs_by_svid, svid):
    """In: the alt allele. Out: PASS, zero error, one full-length forward match."""
    sv = svs_by_svid[svid]
    assert sv.svtype == EXPECTED[svid].svtype

    query = get_reconstructed_query(sv, reference)
    assert len(query) == EXPECTED[svid].reconstructed

    result = assembly_scorer.score_query(query)

    assert result.source is ValidationSource.ASSEMBLY
    assert result.status is QueryValidationStatus.PASS
    assert result.reason is QueryValidationReason.PASS
    assert result.passed is True
    assert result.lowest_error == 0.0
    assert result.best_strand_match == 1
    assert repr(result.cigar) == f'{len(query)}='


@pytest.mark.parametrize('svid', SVIDS, ids=SVID_IDS)
def test_alignment_rejects_reference(assembly_scorer, reference, svs_by_svid, svid):
    """In: the unrearranged reference. Out: FAIL -- the assembly carries the SV instead."""
    query = get_reference_query(svs_by_svid[svid], reference)
    assert len(query) == EXPECTED[svid].reference

    result = assembly_scorer.score_query(query)

    assert result.source is ValidationSource.ASSEMBLY
    assert result.status is QueryValidationStatus.FAIL
    assert result.reason is QueryValidationReason.CIGAR_FAILED  # aligned, then rejected
    assert result.aligned is True
    assert result.passed is False
