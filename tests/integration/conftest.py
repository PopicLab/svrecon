"""Session fixtures: pytest injects these into any test naming one as a parameter, no import."""
from typing import Dict, List

import pytest

from constants import ASSEMBLY_CHROM, ASSEMBLY_FASTA, ASSEMBLY_VCF, CHROM, REFERENCE_FASTA
from generators.generate_small_bam import SimulatedSV, load_sequence, load_simulated_svs


@pytest.fixture(scope='session')
def reference_sequence() -> str:
    return load_sequence(REFERENCE_FASTA, CHROM)


@pytest.fixture(scope='session')
def assembly_sequence() -> str:
    return load_sequence(ASSEMBLY_FASTA, ASSEMBLY_CHROM)


@pytest.fixture(scope='session')
def simulated_svs(reference_sequence, assembly_sequence) -> List[SimulatedSV]:
    return load_simulated_svs(ASSEMBLY_VCF, reference_sequence, assembly_sequence)


@pytest.fixture(scope='session')
def svs_by_svid(simulated_svs) -> Dict[str, SimulatedSV]:
    return {sv.svid: sv for sv in simulated_svs}
