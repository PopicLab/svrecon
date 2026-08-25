"""Session fixtures: pytest injects these into any test naming one as a parameter, no import.

Built once per run and only if a selected test asks -- the FASTA reads, VCF parse and mappy index
happen lazily, not at collection. Dependencies resolve by name: svs_by_svid -> simulated_svs ->
reference, assembly.
"""
import pytest

from generate_small_bam import (ASSEMBLY_CHROM, ASSEMBLY_FASTA, ASSEMBLY_VCF, READS_BAM,
                                load_sequence, load_simulated_svs)
from generate_small_genome import CHROM, REFERENCE_FASTA
from helpers import get_config
from svrecon.scorers.assembly import AssemblyScorer
from svrecon.scorers.reads import ReadScorer


@pytest.fixture(scope='session')
def reference() -> str:
    return load_sequence(REFERENCE_FASTA, CHROM)


@pytest.fixture(scope='session')
def assembly() -> str:
    return load_sequence(ASSEMBLY_FASTA, ASSEMBLY_CHROM)


@pytest.fixture(scope='session')
def simulated_svs(reference, assembly):
    return load_simulated_svs(ASSEMBLY_VCF, reference, assembly)


@pytest.fixture(scope='session')
def svs_by_svid(simulated_svs) -> dict:
    return {sv.svid: sv for sv in simulated_svs}


@pytest.fixture(scope='session')
def read_scorer() -> ReadScorer:
    """ReadScorer over the committed BAM: one supporting read per SV."""
    return ReadScorer(get_config(bam=str(READS_BAM)), chroms={CHROM})


@pytest.fixture(scope='session')
def assembly_scorer(tmp_path_factory) -> AssemblyScorer:
    """AssemblyScorer over the assembly FASTA, with a throwaway mappy index cache."""
    config = get_config(cache_dir=tmp_path_factory.mktemp('mappy_cache'))
    return AssemblyScorer(config, chroms={CHROM}, fasta_path=str(ASSEMBLY_FASTA))
