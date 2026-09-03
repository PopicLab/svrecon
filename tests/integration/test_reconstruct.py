"""Reconstruction against the assembly insilicoSV built: each SV's records in, its allele out.

The records and the reference bytes arrive as conftest fixtures, so the only expected values
here are the ones helpers declares.
"""
from typing import Dict, List

import pytest
from pysam import VariantRecord

from constants import ASSEMBLY_VCF, CHROM, QUERY_BUFFER_SIZE, REFERENCE_FASTA
from helpers import EXPECTED_SVTYPES, SIMULATED_SVIDS, SVID_TEST_IDS
from svrecon.reconstruct import construct_queries
from svrecon.utils import load_fasta_to_bytes, load_grouped_variants_from_vcf


@pytest.fixture(scope='module')
def records_by_svid() -> Dict[str, List[VariantRecord]]:
    return load_grouped_variants_from_vcf(str(ASSEMBLY_VCF))


@pytest.fixture(scope='module')
def reference_bytes() -> Dict[str, bytearray]:
    return load_fasta_to_bytes(str(REFERENCE_FASTA), {CHROM})


@pytest.mark.parametrize('svid', SIMULATED_SVIDS, ids=SVID_TEST_IDS)
def test_reconstruction_matches_assembly(records_by_svid, reference_bytes, assembly_sequence,
                                         svs_by_svid, svid):
    """In: the SV's VCF records. Out: one query, matching the assembly."""
    sv = svs_by_svid[svid]
    assert sv.svtype == EXPECTED_SVTYPES[svid]

    expected_sequence = assembly_sequence[sv.asm_start - QUERY_BUFFER_SIZE:
                                          sv.asm_end + QUERY_BUFFER_SIZE]

    queries = construct_queries(svid, records_by_svid[svid], QUERY_BUFFER_SIZE, reference_bytes)

    assert len(queries) == 1
    assert queries[0].sequence == expected_sequence
