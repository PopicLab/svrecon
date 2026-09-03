"""AssemblyScorer against the assembly insilicoSV built: mappy alignment of each allele."""
import pytest

from constants import ASSEMBLY_FASTA, CHROM
from helpers import (EXPECTED_SVTYPES, SIMULATED_SVIDS, SVID_TEST_IDS, get_true_reconstructed_query,
                     get_reference_query, get_scorer_config)
from svrecon.constants import QueryValidationReason, QueryValidationStatus, ValidationSource
from svrecon.scorers.assembly import AssemblyScorer


@pytest.fixture(scope='module')
def assembly_scorer(tmp_path_factory) -> AssemblyScorer:
    """AssemblyScorer over the assembly FASTA, with a throwaway mappy index cache."""
    config = get_scorer_config(cache_dir=tmp_path_factory.mktemp('mappy_cache'))
    return AssemblyScorer(config, chroms={CHROM}, fasta_path=str(ASSEMBLY_FASTA))


@pytest.mark.parametrize('svid', SIMULATED_SVIDS, ids=SVID_TEST_IDS)
def test_alignment_supports_reconstructed(assembly_scorer, assembly_sequence,
                                          svs_by_svid, svid):
    """In: the alt allele. Out: PASS, zero error, one full-length forward match."""
    sv = svs_by_svid[svid]
    assert sv.svtype == EXPECTED_SVTYPES[svid]

    query = get_true_reconstructed_query(sv, assembly_sequence)

    result = assembly_scorer.score_query(query)

    assert result.source is ValidationSource.ASSEMBLY
    assert result.status is QueryValidationStatus.PASS
    assert result.reason is QueryValidationReason.PASS
    assert result.passed is True
    assert result.lowest_error == 0.0
    assert result.best_strand_match == 1
    assert repr(result.cigar) == f'{len(query)}='


@pytest.mark.parametrize('svid', SIMULATED_SVIDS, ids=SVID_TEST_IDS)
def test_alignment_rejects_reference(assembly_scorer, reference_sequence,
                                     svs_by_svid, svid):
    """In: the unrearranged reference. Out: FAIL -- the assembly carries the SV instead."""
    query = get_reference_query(svs_by_svid[svid], reference_sequence)

    result = assembly_scorer.score_query(query)

    assert result.source is ValidationSource.ASSEMBLY
    assert result.status is QueryValidationStatus.FAIL
    assert result.reason is QueryValidationReason.CIGAR_FAILED  # aligned, then rejected
    assert result.aligned is True
    assert result.passed is False
