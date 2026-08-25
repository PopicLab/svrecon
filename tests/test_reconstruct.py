"""Reconstruction against the assembly insilicoSV built: the rebuilt allele must match it exactly."""
import pytest

from generate_small_bam import ASSEMBLY_VCF
from generate_small_genome import CHROM, REFERENCE_FASTA
from helpers import BUFFER, EXPECTED, SVID_IDS, SVIDS
from svrecon.reconstruct import simulate_subsequences
from svrecon.utils import group_variants_by_id, load_fasta_to_bytes


@pytest.fixture(scope='module')
def records_by_svid() -> dict:
    return group_variants_by_id(str(ASSEMBLY_VCF))


@pytest.fixture(scope='module')
def reference_bytes() -> dict:
    return load_fasta_to_bytes(str(REFERENCE_FASTA), {CHROM})


def get_expected_sequence(sv, reference: str) -> str:
    """The allele insilicoSV built, with BUFFER bp of reference each side."""
    return (reference[sv.ref_start - BUFFER:sv.ref_start] + sv.alt_sequence
            + reference[sv.ref_end:sv.ref_end + BUFFER])


@pytest.mark.parametrize('svid', SVIDS, ids=SVID_IDS)
def test_reconstruction_matches_assembly(records_by_svid, reference_bytes, reference,
                                         svs_by_svid, svid):
    """In: the SV's VCF records. Out: one window, matching the assembly."""
    sv = svs_by_svid[svid]
    assert sv.svtype == EXPECTED[svid].svtype

    queries = simulate_subsequences(records_by_svid[svid], BUFFER, reference_bytes)

    assert len(queries) == 1
    assert queries[0].sequence == get_expected_sequence(sv, reference)
