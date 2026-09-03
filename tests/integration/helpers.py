"""Expected values for the simulated fixture, and the query builders the scorer tests align."""
from types import SimpleNamespace

from constants import CHROM, PLACEHOLDER_GRAMMAR, QUERY_BUFFER_SIZE
from svrecon.config import DEFAULTS
from svrecon.reconstruct import Query


# Keyed by SVID
EXPECTED_SVTYPES = {
    'sv0':  'DEL',
    'sv1':  'DEL',
    'sv2':  'INV',
    'sv3':  'INV',
    'sv4':  'DUP',
    'sv5':  'DUP',
    'sv6':  'DUP_INV',
    'sv7':  'DUP_INV',
    'sv8':  'delINV',
    'sv9':  'delINV',
    'sv10': 'INVdel',
    'sv11': 'INVdel',
    'sv12': 'delINVdel',
    'sv13': 'delINVdel',
    'sv14': 'delINVdup',
    'sv15': 'delINVdup',
    'sv16': 'dupINVdel',
    'sv17': 'dupINVdel',
    'sv18': 'dupINVdup',
    'sv19': 'dupINVdup',
}
SIMULATED_SVIDS = list(EXPECTED_SVTYPES)
SVID_TEST_IDS = [f'{svid}-{svtype}' for svid, svtype in EXPECTED_SVTYPES.items()]


def get_scorer_config(**overrides) -> SimpleNamespace:
    """Stand in config, with standard config defaults"""
    return SimpleNamespace(**{**DEFAULTS, **overrides})

def get_true_reconstructed_query(sv, assembly_sequence: str) -> Query:
    """Query of sv's true assembly sequence"""
    query_ref_start, query_ref_end = sv.ref_start - QUERY_BUFFER_SIZE, sv.ref_end + QUERY_BUFFER_SIZE
    return Query(chrom=CHROM, svtype=sv.svtype, svid=sv.svid, grammar=PLACEHOLDER_GRAMMAR,
                 sequence=assembly_sequence[sv.asm_start - QUERY_BUFFER_SIZE:
                                            sv.asm_end + QUERY_BUFFER_SIZE],
                 ref_start=query_ref_start, ref_end=query_ref_end,
                 recon_segments=[], ref_segments=[(query_ref_start, query_ref_end)], ref_sequence='',
                 buffer=QUERY_BUFFER_SIZE)


def get_reference_query(sv, reference_sequence: str) -> Query:
    """Query constructed with ref sequence where sv was derived from"""
    query_ref_start, query_ref_end = sv.ref_start - QUERY_BUFFER_SIZE, sv.ref_end + QUERY_BUFFER_SIZE
    return Query(chrom=CHROM, svtype=sv.svtype, svid=sv.svid, grammar=PLACEHOLDER_GRAMMAR,
                 sequence=reference_sequence[query_ref_start:query_ref_end],
                 ref_start=query_ref_start, ref_end=query_ref_end,
                 recon_segments=[], ref_segments=[(query_ref_start, query_ref_end)], ref_sequence='',
                 buffer=QUERY_BUFFER_SIZE)
