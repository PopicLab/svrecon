"""Reconstruction tests."""
from typing import List, NamedTuple

import pytest

from generate_small_bam import ASSEMBLY_VCF
from generate_small_genome import CHROM, REFERENCE_FASTA
from helpers import BUFFER as ASSEMBLY_BUFFER  # bp count, not the flank sequence BUFFER below
from helpers import EXPECTED, SVID_IDS, SVIDS
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

    queries = construct_queries(records_by_svid[svid], ASSEMBLY_BUFFER, reference_bytes)

    assert len(queries) == 1
    assert queries[0].sequence == get_expected_sequence(sv, reference, ASSEMBLY_BUFFER)


# --- hand-derived from grammar, not insilicoSV output -------------------------------------------

SYN_CHROM = 'chrT'
BUFFER = 'GT'
BUFFER_SIZE = len(BUFFER)  # every grammar case buffers by exactly the flank it carries each side
A, B, C = 'AAAA', 'CACACA', 'ACACACAC'                                         # A/C only
a, b, c = reverse_complement(A), reverse_complement(B), reverse_complement(C)  # T/G only

SHORT_DISP = 'CCAA' # exactly half of buffer size
LONG_DISP = SHORT_DISP * 2
DISP_HEAD, DISP_TAIL = LONG_DISP[:BUFFER_SIZE], LONG_DISP[-BUFFER_SIZE:]


class DummyRecord:
    """Duck-typed stand-in for pysam.VariantRecord: only the fields reconstruct.py reads."""
    def __init__(self, start, stop, alt, target=None, insord=-1):
        self.chrom = SYN_CHROM
        self.start = start
        self.stop = stop
        self.alts = (f'<{alt}>',)
        self.info = {'SVID': 'sv', 'SVTYPE': 'sv'}
        if target is not None:
            self.info['TARGET'], self.info['INSORD'] = target, insord


class GrammarCase(NamedTuple):
    name: str
    reference: str                 #  BUFFER + contig + BUFFER
    buffer: int
    records: List[DummyRecord]
    expected_sequences: List[str]  # [buffer + contig + buffer]

GRAMMAR_CASES = [
    GrammarCase('DEL', BUFFER + A + BUFFER, BUFFER_SIZE, [
        DummyRecord(start=len(BUFFER), stop=len(BUFFER + A), alt='DEL')], # A -> 
        [BUFFER + BUFFER]),
    GrammarCase('INV', BUFFER + A + BUFFER, BUFFER_SIZE, [
        DummyRecord(start=len(BUFFER), stop=len(BUFFER + A), alt='INV')],# A -> a
        [BUFFER + a + BUFFER]),
    GrammarCase('DUP', BUFFER + A + BUFFER, BUFFER_SIZE, [
        DummyRecord(start=len(BUFFER), stop=len(BUFFER + A), alt='DUP')], # A -> AA
        [BUFFER + A + A + BUFFER]),
    GrammarCase('INV_DUP', BUFFER + A + BUFFER, BUFFER_SIZE, [
        DummyRecord(start=len(BUFFER), stop=len(BUFFER + A), alt='INV_DUP', target=len(BUFFER + A), insord=0)], # A -> Aa
        [BUFFER + A + a + BUFFER]),
    GrammarCase('DUP_INV', BUFFER + A + BUFFER, BUFFER_SIZE, [
        DummyRecord(start=len(BUFFER), stop=len(BUFFER + A), alt='INV'), #A -> a
        DummyRecord(start=len(BUFFER), stop=len(BUFFER + A), alt='COPYinv-PASTE', target=len(BUFFER + A), insord=0)], # A -> aa
        [BUFFER + a + a + BUFFER]),
    GrammarCase('delINV', BUFFER + A + B + BUFFER, BUFFER_SIZE, [
        DummyRecord(start=len(BUFFER + A), stop=len(BUFFER + A + B), alt='INV'), # AB -> Ab
        DummyRecord(start=len(BUFFER), stop=len(BUFFER + A), alt='CUT')], # AB -> b
        [BUFFER + b + BUFFER]),
    GrammarCase('INVdel', BUFFER + A + B + BUFFER, BUFFER_SIZE, [
        DummyRecord(start=len(BUFFER + A), stop=len(BUFFER + A + B), alt='CUT'), # AB -> A
        DummyRecord(start=len(BUFFER), stop=len(BUFFER + A), alt='INV')], # AB -> a
        [BUFFER + a + BUFFER]),
    GrammarCase('dupINV', BUFFER + A + B + BUFFER, BUFFER_SIZE, [
        DummyRecord(start=len(BUFFER + A), stop=len(BUFFER + A + B), alt='INV'), # AB -> Ab
        DummyRecord(start=len(BUFFER), stop=len(BUFFER + A), alt='COPYinv-PASTE', target=len(BUFFER + A + B), insord=0)], # AB -> Aba
        [BUFFER + A + b + a + BUFFER]),
    GrammarCase('INVdup', BUFFER + A + B + BUFFER, BUFFER_SIZE, [
        DummyRecord(start=len(BUFFER + A), stop=len(BUFFER + A + B), alt='COPY-PASTE', target=len(BUFFER + A + B), insord=1), # AB -> ABB
        DummyRecord(start=len(BUFFER), stop=len(BUFFER + A), alt='CUTinv-PASTE', target=len(BUFFER + A + B), insord=0), # invPaste AB -> ABaB, CUT AB-> AbaB, INV AB-> baB
        DummyRecord(start=len(BUFFER + A), stop=len(BUFFER + A + B), alt='INV')], 
        [BUFFER + b + a + B + BUFFER]),
    GrammarCase('dupINVdup', BUFFER + A + B + C + BUFFER, BUFFER_SIZE, [
        DummyRecord(start=len(BUFFER + A + B), stop=len(BUFFER + A + B + C), alt='INV'), # ABC -> ABc
        DummyRecord(start=len(BUFFER + A + B), stop=len(BUFFER + A + B + C), alt='COPY-PASTE', target=len(BUFFER + A + B + C), insord=2), # ABC -> ABcC
        DummyRecord(start=len(BUFFER), stop=len(BUFFER + A), alt='COPYinv-PASTE', target=len(BUFFER + A + B + C), insord=1), # ABC -> ABcaC
        DummyRecord(start=len(BUFFER + A), stop=len(BUFFER + A + B), alt='CUTinv-PASTE', target=len(BUFFER + A + B + C), insord=0)], # ABC -> AcbaC
        [BUFFER + A + c + b + a + C + BUFFER]),
    GrammarCase('delINVdel', BUFFER + A + B + C + BUFFER, BUFFER_SIZE, [
        DummyRecord(start=len(BUFFER + A + B), stop=len(BUFFER + A + B + C), alt='CUT'), # ABC -> AB
        DummyRecord(start=len(BUFFER + A), stop=len(BUFFER + A + B), alt='INV'), # ABC -> Ab
        DummyRecord(start=len(BUFFER), stop=len(BUFFER + A), alt='CUT')], # ABC -> b
        [BUFFER + b + BUFFER]),
    GrammarCase('delINVdup', BUFFER + A + B + C + BUFFER, BUFFER_SIZE, [
        DummyRecord(start=len(BUFFER + A + B), stop=len(BUFFER + A + B + C), alt='COPY-PASTE', target=len(BUFFER + A + B + C), insord=1), # ABC -> ABCC
        DummyRecord(start=len(BUFFER + A), stop=len(BUFFER + A + B), alt='CUTinv-PASTE', target=len(BUFFER + A + B + C), insord=0), # ABC-> ABCbC ->ABcbC -> AcbC
        DummyRecord(start=len(BUFFER + A + B), stop=len(BUFFER + A + B + C), alt='INV'), 
        DummyRecord(start=len(BUFFER), stop=len(BUFFER + A), alt='CUT')], # ABC -> cbC
        [BUFFER + c + b + C + BUFFER]),
    GrammarCase('dupINVdel', BUFFER + A + B + C + BUFFER, BUFFER_SIZE, [
        DummyRecord(start=len(BUFFER + A + B), stop=len(BUFFER + A + B + C), alt='CUT'), # ABC -> AB
        DummyRecord(start=len(BUFFER), stop=len(BUFFER + A), alt='COPYinv-PASTE', target=len(BUFFER + A + B), insord=0), # ABC -> ABa
        DummyRecord(start=len(BUFFER + A), stop=len(BUFFER + A + B), alt='INV')], # ABC -> Aba
        [BUFFER + A + b + a + BUFFER]),

    # dispersed, twice each: SHORT_DISP keeps one window, LONG_DISP splits into two
    GrammarCase('dDUP-short', BUFFER + A + SHORT_DISP + BUFFER, BUFFER_SIZE, [
        DummyRecord(start=len(BUFFER), stop=len(BUFFER + A), alt='dDUP', target=len(BUFFER + A + SHORT_DISP), insord=0)],
        [BUFFER + A + SHORT_DISP + A + BUFFER]),
    GrammarCase('dDUP-long', BUFFER + A + LONG_DISP + BUFFER, BUFFER_SIZE, [
        DummyRecord(start=len(BUFFER), stop=len(BUFFER + A), alt='dDUP', target=len(BUFFER + A + LONG_DISP), insord=0)],
        [BUFFER + A + DISP_HEAD,
         DISP_TAIL + A + BUFFER]),
    GrammarCase('INV_dDUP-short', BUFFER + A + SHORT_DISP + BUFFER, BUFFER_SIZE, [
        DummyRecord(start=len(BUFFER), stop=len(BUFFER + A), alt='INV_dDUP', target=len(BUFFER + A + SHORT_DISP), insord=0)],
        [BUFFER + A + SHORT_DISP + a + BUFFER]),
    GrammarCase('INV_dDUP-long', BUFFER + A + LONG_DISP + BUFFER, BUFFER_SIZE, [
        DummyRecord(start=len(BUFFER), stop=len(BUFFER + A), alt='INV_dDUP', target=len(BUFFER + A + LONG_DISP), insord=0)],
        [BUFFER + A + DISP_HEAD,
         DISP_TAIL + a + BUFFER]),
    GrammarCase('dDUP_INV-short', BUFFER + A + SHORT_DISP + BUFFER, BUFFER_SIZE, [
        DummyRecord(start=len(BUFFER), stop=len(BUFFER + A), alt='INV'),
        DummyRecord(start=len(BUFFER), stop=len(BUFFER + A), alt='COPYinv-PASTE', target=len(BUFFER + A + SHORT_DISP), insord=0)],
        [BUFFER + a + SHORT_DISP + a + BUFFER]),
    GrammarCase('dDUP_INV-long', BUFFER + A + LONG_DISP + BUFFER, BUFFER_SIZE, [
        DummyRecord(start=len(BUFFER), stop=len(BUFFER + A), alt='INV'),
        DummyRecord(start=len(BUFFER), stop=len(BUFFER + A), alt='COPYinv-PASTE', target=len(BUFFER + A + LONG_DISP), insord=0)],
        [BUFFER + a + DISP_HEAD,
         DISP_TAIL + a + BUFFER]),
    GrammarCase('INV_nrTRA-short', BUFFER + A + SHORT_DISP + BUFFER, BUFFER_SIZE, [
        DummyRecord(start=len(BUFFER), stop=len(BUFFER + A), alt='INV_nrTRA', target=len(BUFFER + A + SHORT_DISP), insord=0)],
        [BUFFER + SHORT_DISP + a + BUFFER]),
    GrammarCase('INV_nrTRA-long', BUFFER + A + LONG_DISP + BUFFER, BUFFER_SIZE, [
        DummyRecord(start=len(BUFFER), stop=len(BUFFER + A), alt='INV_nrTRA', target=len(BUFFER + A + LONG_DISP), insord=0)],
        [BUFFER + DISP_HEAD,
         DISP_TAIL + a + BUFFER]),
    GrammarCase('nrTRA-short', BUFFER + A + SHORT_DISP + BUFFER, BUFFER_SIZE, [
        DummyRecord(start=len(BUFFER), stop=len(BUFFER + A), alt='nrTRA', target=len(BUFFER + A + SHORT_DISP), insord=0)],
        [BUFFER + SHORT_DISP + A + BUFFER]),
    GrammarCase('nrTRA-long', BUFFER + A + LONG_DISP + BUFFER, BUFFER_SIZE, [
        DummyRecord(start=len(BUFFER), stop=len(BUFFER + A), alt='nrTRA', target=len(BUFFER + A + LONG_DISP), insord=0)],
        [BUFFER + DISP_HEAD,
         DISP_TAIL + A + BUFFER]),
    GrammarCase('rTRA-short', BUFFER + A + SHORT_DISP + B + BUFFER, BUFFER_SIZE, [
        DummyRecord(start=len(BUFFER), stop=len(BUFFER + A), alt='CUT-PASTE', target=len(BUFFER + A + SHORT_DISP), insord=1), # 
        DummyRecord(start=len(BUFFER + A + SHORT_DISP), stop=len(BUFFER + A + SHORT_DISP + B), alt='CUT-PASTE', target=len(BUFFER), insord=0)],
        [BUFFER + B + SHORT_DISP + A + BUFFER]),
    GrammarCase('rTRA-long', BUFFER + A + LONG_DISP + B + BUFFER, BUFFER_SIZE, [
        DummyRecord(start=len(BUFFER), stop=len(BUFFER + A), alt='CUT-PASTE', target=len(BUFFER + A + LONG_DISP), insord=1),
        DummyRecord(start=len(BUFFER + A + LONG_DISP), stop=len(BUFFER + A + LONG_DISP + B), alt='CUT-PASTE', target=len(BUFFER), insord=0)],
        [BUFFER + B + DISP_HEAD,
         DISP_TAIL + A + BUFFER]),
    GrammarCase('INV_rTRA-short', BUFFER + A + SHORT_DISP + B + BUFFER, BUFFER_SIZE, [
        DummyRecord(start=len(BUFFER), stop=len(BUFFER + A), alt='CUTinv-PASTE', target=len(BUFFER + A + SHORT_DISP), insord=1),
        DummyRecord(start=len(BUFFER + A + SHORT_DISP), stop=len(BUFFER + A + SHORT_DISP + B), alt='CUTinv-PASTE', target=len(BUFFER), insord=0)],
        [BUFFER + b + SHORT_DISP + a + BUFFER]),
    GrammarCase('INV_rTRA-long', BUFFER + A + LONG_DISP + B + BUFFER, BUFFER_SIZE, [
        DummyRecord(start=len(BUFFER), stop=len(BUFFER + A), alt='CUTinv-PASTE', target=len(BUFFER + A + LONG_DISP), insord=1),
        DummyRecord(start=len(BUFFER + A + LONG_DISP), stop=len(BUFFER + A + LONG_DISP + B), alt='CUTinv-PASTE', target=len(BUFFER), insord=0)],
        [BUFFER + b + DISP_HEAD,
         DISP_TAIL + a + BUFFER]),
    GrammarCase('INS_iDEL-short', BUFFER + A + SHORT_DISP + B + BUFFER, BUFFER_SIZE, [
        DummyRecord(start=len(BUFFER + A + SHORT_DISP), stop=len(BUFFER + A + SHORT_DISP + B), alt='CUT'),
        DummyRecord(start=len(BUFFER), stop=len(BUFFER + A), alt='CUT-PASTE', target=len(BUFFER + A + SHORT_DISP), insord=0)],
        [BUFFER + SHORT_DISP + A + BUFFER]),
    GrammarCase('INS_iDEL-long', BUFFER + A + LONG_DISP + B + BUFFER, BUFFER_SIZE, [
        DummyRecord(start=len(BUFFER), stop=len(BUFFER + A), alt='CUT-PASTE', target=len(BUFFER + A + LONG_DISP), insord=0),
        DummyRecord(start=len(BUFFER + A + LONG_DISP), stop=len(BUFFER + A + LONG_DISP + B), alt='CUT')],
        [BUFFER + DISP_HEAD,
         DISP_TAIL + A + BUFFER]),
    GrammarCase('dDUP_iDEL-short', BUFFER + A + SHORT_DISP + B + BUFFER, BUFFER_SIZE, [
        DummyRecord(start=len(BUFFER), stop=len(BUFFER + A), alt='COPY-PASTE', target=len(BUFFER + A + SHORT_DISP), insord=0),
        DummyRecord(start=len(BUFFER + A + SHORT_DISP), stop=len(BUFFER + A + SHORT_DISP + B), alt='CUT')],
        [BUFFER + A + SHORT_DISP + A + BUFFER]),
    GrammarCase('dDUP_iDEL-long', BUFFER + A + LONG_DISP + B + BUFFER, BUFFER_SIZE, [
        DummyRecord(start=len(BUFFER + A + LONG_DISP), stop=len(BUFFER + A + LONG_DISP + B), alt='CUT'),
        DummyRecord(start=len(BUFFER), stop=len(BUFFER + A), alt='COPY-PASTE', target=len(BUFFER + A + LONG_DISP), insord=0)],
        [BUFFER + A + DISP_HEAD,
         DISP_TAIL + A + BUFFER]),
]


@pytest.mark.parametrize('reference,buffer,records,expected_sequences',
                         [case[1:] for case in GRAMMAR_CASES],
                         ids=[case.name for case in GRAMMAR_CASES])
def test_reconstruction_matches_grammar(reference, buffer, records, expected_sequences):
    """Records built directly from each type's known grammar."""
    ref = {SYN_CHROM: reference.encode()}

    queries = construct_queries(records, buffer, ref)

    assert [query.sequence for query in queries] == expected_sequences
