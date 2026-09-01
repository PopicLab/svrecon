"""Reconstruction against the assembly insilicoSV built: the rebuilt allele must match it exactly."""
import pytest

from generate_small_bam import ASSEMBLY_VCF
from generate_small_genome import CHROM, REFERENCE_FASTA
from helpers import BUFFER, EXPECTED, SVID_IDS, SVIDS
from svrecon.reconstruct import construct_queries
from svrecon.utils import group_variants_by_id, load_fasta_to_bytes, reverse_complement


@pytest.fixture(scope='module')
def records_by_svid() -> dict:
    return group_variants_by_id(str(ASSEMBLY_VCF))


@pytest.fixture(scope='module')
def reference_bytes() -> dict:
    return load_fasta_to_bytes(str(REFERENCE_FASTA), {CHROM})


def get_expected_sequence(sv, reference: str, buffer: int) -> str:
    """The allele insilicoSV built, with ``buffer`` bp of reference each side."""
    return (reference[sv.ref_start - buffer:sv.ref_start] + sv.alt_sequence
            + reference[sv.ref_end:sv.ref_end + buffer])


@pytest.mark.parametrize('svid', SVIDS, ids=SVID_IDS)
def test_reconstruction_matches_assembly(records_by_svid, reference_bytes, reference,
                                         svs_by_svid, svid):
    """In: the SV's VCF records. Out: one window, matching the assembly."""
    sv = svs_by_svid[svid]
    assert sv.svtype == EXPECTED[svid].svtype

    queries = construct_queries(records_by_svid[svid], BUFFER, reference_bytes)

    assert len(queries) == 1
    assert queries[0].sequence == get_expected_sequence(sv, reference, BUFFER)


# --- hand-derived from grammar, not insilicoSV output -------------------------------------------

SYN_CHROM, FLANK = 'chrT', 'GT'
A, B, C = 'AAAA', 'CACACA', 'ACACACAC'                                   # letters A/C only
a, b, c = reverse_complement(A), reverse_complement(B), reverse_complement(C)  # letters T/G only


class _DummyRecord:
    """Duck-typed stand-in for pysam.VariantRecord: only the fields reconstruct.py reads."""
    def __init__(self, start, stop, alt, target=None, insord=-1):
        self.chrom = SYN_CHROM
        self.start = start
        self.stop = stop
        self.alts = (f'<{alt}>',)
        self.info = {'SVID': 'sv', 'SVTYPE': 'sv'}
        if target is not None:
            self.info['TARGET'], self.info['INSORD'] = target, insord


# reference (A/B/C back to back), ops (alt, start, stop, target, insord), expected core
GRAMMAR_CASES = [
    ('DEL', A, [
        ('DEL', 0, 4, None, None)],
     ''),
    ('INV', A, [
        ('INV', 0, 4, None, None)],
     a),
    ('DUP', A, [
        ('DUP', 0, 4, None, None)],
     A + A),
    ('DUP_INV', A, [
        ('INV', 0, 4, None, None),
        ('COPYinv-PASTE', 0, 4, 4, 0)],
     a + a),
    ('delINV', A + B, [
        ('CUT', 0, 4, None, None),
        ('INV', 4, 10, None, None)],
     b),
    ('INVdel', A + B, [
        ('INV', 0, 4, None, None),
        ('CUT', 4, 10, None, None)],
     a),
    ('delINVdel', A + B + C, [
        ('CUT', 0, 4, None, None),
        ('INV', 4, 10, None, None),
        ('CUT', 10, 18, None, None)],
     b),
    ('delINVdup', A + B + C, [
        ('CUT', 0, 4, None, None),
        ('CUTinv-PASTE', 4, 10, 18, 0),
        ('INV', 10, 18, None, None),
        ('COPY-PASTE', 10, 18, 18, 1)],
     c + b + C),
    ('dupINVdel', A + B + C, [
        ('COPYinv-PASTE', 0, 4, 10, 0),
        ('INV', 4, 10, None, None),
        ('CUT', 10, 18, None, None)],
     A + b + a),
    ('dupINVdup', A + B + C, [
        ('COPYinv-PASTE', 0, 4, 18, 1),
        ('CUTinv-PASTE', 4, 10, 18, 0),
        ('INV', 10, 18, None, None),
        ('COPY-PASTE', 10, 18, 18, 2)],
     A + c + b + a + C),
]


@pytest.mark.parametrize('core_ref,ops,expected_core', [case[1:] for case in GRAMMAR_CASES],
                         ids=[case[0] for case in GRAMMAR_CASES])
def test_reconstruction_matches_grammar(core_ref, ops, expected_core):
    """Records built directly from each type's known grammar, not insilicoSV's own output."""
    buffer = len(FLANK)
    ref = {SYN_CHROM: (FLANK + core_ref + FLANK).encode()}
    records = [_DummyRecord(start + buffer, stop + buffer, alt,
                            None if target is None else target + buffer, insord)
              for alt, start, stop, target, insord in ops]

    queries = construct_queries(records, buffer, ref)

    assert len(queries) == 1
    assert queries[0].sequence == FLANK + expected_core + FLANK
