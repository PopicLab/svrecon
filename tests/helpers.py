"""Expected values and query builders"""
from types import SimpleNamespace
from typing import NamedTuple

from generate_small_genome import CHROM
from svrecon.config import DEFAULTS
from svrecon.reconstruct import Query

BUFFER = 500  # reference context each side of the reconstruction; under the fixture's READ_BUFFER


class ExpectedSV(NamedTuple):
    """One simulated SV: its type, and each allele's query length (2 * BUFFER + haplotype)."""
    svtype: str
    reconstructed: int
    reference: int


# Keyed by SVID, which insilicoSV assigns in variant_sets order -- two per type, so sv0/sv1 are
# the DELs, sv2/sv3 the INVs, and so on. Declared, not discovered: parametrize needs them at
# collection, and pinning svtype here catches the mapping shifting under us.
EXPECTED = {
    'sv0':  ExpectedSV('DEL',       1_000, 1_300),
    'sv1':  ExpectedSV('DEL',       1_000, 1_300),
    'sv2':  ExpectedSV('INV',       1_300, 1_300),
    'sv3':  ExpectedSV('INV',       1_300, 1_300),
    'sv4':  ExpectedSV('DUP',       1_600, 1_300),
    'sv5':  ExpectedSV('DUP',       1_600, 1_300),
    'sv6':  ExpectedSV('DUP_INV',   1_600, 1_300),
    'sv7':  ExpectedSV('DUP_INV',   1_600, 1_300),
    'sv8':  ExpectedSV('delINV',    1_500, 1_800),
    'sv9':  ExpectedSV('delINV',    1_500, 1_800),
    'sv10': ExpectedSV('INVdel',    1_300, 1_800),
    'sv11': ExpectedSV('INVdel',    1_300, 1_800),
    'sv12': ExpectedSV('delINVdel', 1_500, 2_200),
    'sv13': ExpectedSV('delINVdel', 1_500, 2_200),
    'sv14': ExpectedSV('delINVdup', 2_300, 2_200),
    'sv15': ExpectedSV('delINVdup', 2_300, 2_200),
    'sv16': ExpectedSV('dupINVdel', 2_100, 2_200),
    'sv17': ExpectedSV('dupINVdel', 2_100, 2_200),
    'sv18': ExpectedSV('dupINVdup', 2_900, 2_200),
    'sv19': ExpectedSV('dupINVdup', 2_900, 2_200),
}
SVIDS = list(EXPECTED)
SVID_IDS = [f'{svid}-{expected.svtype}' for svid, expected in EXPECTED.items()]  # readable test ids


def get_config(**overrides) -> SimpleNamespace:
    """DEFAULTS plus overrides. Avoids Config, whose init reconfigures logging."""
    return SimpleNamespace(**{**DEFAULTS, **overrides})


GRAMMAR = 'dummy'  # scorers never read it; only construct_queries derives a real one


def get_reconstructed_query(sv, reference: str) -> Query:
    """Query from the simulated SV -- positive testing."""
    ref_start, ref_end = sv.ref_start - BUFFER, sv.ref_end + BUFFER
    sequence = (reference[ref_start:sv.ref_start] + sv.alt_sequence
                + reference[sv.ref_end:ref_end])
    alt_end = BUFFER + len(sv.alt_sequence)
    return Query(chrom=CHROM, svtype=sv.svtype, svid=sv.svid, grammar=GRAMMAR, sequence=sequence,
                 ref_start=ref_start, ref_end=ref_end,
                 recon_segments=[(0, BUFFER), (BUFFER, alt_end), (alt_end, len(sequence))],
                 ref_segments=[(ref_start, sv.ref_start), (sv.ref_start, sv.ref_end),
                               (sv.ref_end, ref_end)],
                 ref_sequence=reference[ref_start:ref_end], buffer=BUFFER)


def get_reference_query(sv, reference: str) -> Query:
    """The untouched reference over the same span -- negative testing."""
    ref_start, ref_end = sv.ref_start - BUFFER, sv.ref_end + BUFFER
    sequence = reference[ref_start:ref_end]
    return Query(chrom=CHROM, svtype=sv.svtype, svid=sv.svid, grammar=GRAMMAR, sequence=sequence,
                 ref_start=ref_start, ref_end=ref_end,
                 recon_segments=[(0, len(sequence))], ref_segments=[(ref_start, ref_end)],
                 ref_sequence=sequence, buffer=BUFFER)
