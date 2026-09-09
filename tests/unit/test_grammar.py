"""Reconstruction against each SV type's known grammar -- hand-derived, not insilicoSV output.

Needs no fixture data: every case declares its own tiny reference, its VCF records, and the
sequences reconstruction must produce from them.
"""
from typing import List, NamedTuple

import pytest

from svrecon.reconstruct import construct_queries
from svrecon.utils import group_records_by_id, reverse_complement


SYNTHETIC_CHROM = 'chrT'
BUFFER = 'GT'
BUFFER_SIZE = len(BUFFER)  # every grammar case buffers by exactly the flank it carries each side
A, B, C = 'AAAA', 'CACACA', 'ACACACAC'                                         # A/C only
a, b, c = reverse_complement(A), reverse_complement(B), reverse_complement(C)  # T/G only

SHORT_DISP = 'CCAA' # exactly half of buffer size
LONG_DISP = SHORT_DISP * 2
DISP_HEAD, DISP_TAIL = LONG_DISP[:BUFFER_SIZE], LONG_DISP[-BUFFER_SIZE:]


class DummyRecord:
    """realistic input of a stripped down variant record."""
    def __init__(self, start, stop, sv_type=None, op_type=None, svid=None, target=None, insord=-1):
        self.chrom = SYNTHETIC_CHROM
        self.start = start
        self.stop = stop
        self.info = {'SVTYPE': sv_type}
        if op_type:
            self.info['OP_TYPE'] = op_type
        if svid is not None:
            self.info['SVID'] = 'dummyid'
        if target is not None:
            self.info['TARGET'], self.info['INSORD'] = target, insord


class GrammarCase(NamedTuple):
    name: str
    reference: str                 #  BUFFER + contig + BUFFER
    buffer_size: int
    records: List[DummyRecord]
    expected_sequences: List[str]  # [[buffer + contig + buffer],...]
GRAMMAR_CASES = [
    GrammarCase('DEL', BUFFER + A + BUFFER, BUFFER_SIZE, [
        DummyRecord(start=len(BUFFER), stop=len(BUFFER + A), sv_type='DEL')], # A ->
        [BUFFER + BUFFER]),
    GrammarCase('INV', BUFFER + A + BUFFER, BUFFER_SIZE, [
        DummyRecord(start=len(BUFFER), stop=len(BUFFER + A), sv_type='INV')],# A -> a
        [BUFFER + a + BUFFER]),
    GrammarCase('DUP', BUFFER + A + BUFFER, BUFFER_SIZE, [
        DummyRecord(start=len(BUFFER), stop=len(BUFFER + A), sv_type='DUP')], # A -> AA
        [BUFFER + A + A + BUFFER]),
    GrammarCase('INV_DUP', BUFFER + A + BUFFER, BUFFER_SIZE, [
        DummyRecord(start=len(BUFFER), stop=len(BUFFER + A), sv_type='INV_DUP', target=len(BUFFER + A), insord=0)], # A -> Aa
        [BUFFER + A + a + BUFFER]),
    GrammarCase('DUP_INV', BUFFER + A + BUFFER, BUFFER_SIZE, [
        DummyRecord(start=len(BUFFER), stop=len(BUFFER + A), op_type='INV', svid='DUP_INV'), #A -> a
        DummyRecord(start=len(BUFFER), stop=len(BUFFER + A), op_type='COPYinv-PASTE', svid='DUP_INV', target=len(BUFFER + A), insord=0)], # A -> aa
        [BUFFER + a + a + BUFFER]),
    GrammarCase('delINV', BUFFER + A + B + BUFFER, BUFFER_SIZE, [
        DummyRecord(start=len(BUFFER + A), stop=len(BUFFER + A + B), op_type='INV', svid='delINV'), # AB -> Ab
        DummyRecord(start=len(BUFFER), stop=len(BUFFER + A), op_type='CUT', svid='delINV')], # AB -> b
        [BUFFER + b + BUFFER]),
    GrammarCase('INVdel', BUFFER + A + B + BUFFER, BUFFER_SIZE, [
        DummyRecord(start=len(BUFFER + A), stop=len(BUFFER + A + B), op_type='CUT', svid='INVdel'), # AB -> A
        DummyRecord(start=len(BUFFER), stop=len(BUFFER + A), op_type='INV', svid='INVdel')], # AB -> a
        [BUFFER + a + BUFFER]),
    GrammarCase('dupINV', BUFFER + A + B + BUFFER, BUFFER_SIZE, [
        DummyRecord(start=len(BUFFER + A), stop=len(BUFFER + A + B), op_type='INV', svid='dupINV'), # AB -> Ab
        DummyRecord(start=len(BUFFER), stop=len(BUFFER + A), op_type='COPYinv-PASTE', svid='dupINV', target=len(BUFFER + A + B), insord=0)], # AB -> Aba
        [BUFFER + A + b + a + BUFFER]),
    GrammarCase('INVdup', BUFFER + A + B + BUFFER, BUFFER_SIZE, [
        DummyRecord(start=len(BUFFER + A), stop=len(BUFFER + A + B), op_type='COPY-PASTE', svid='INVdup', target=len(BUFFER + A + B), insord=1), # AB -> ABB
        DummyRecord(start=len(BUFFER), stop=len(BUFFER + A), op_type='CUTinv-PASTE', svid='INVdup', target=len(BUFFER + A + B), insord=0), # invPaste AB -> ABaB, CUT AB-> AbaB, INV AB-> baB
        DummyRecord(start=len(BUFFER + A), stop=len(BUFFER + A + B), op_type='INV', svid='INVdup')], 
        [BUFFER + b + a + B + BUFFER]),
    GrammarCase('dupINVdup', BUFFER + A + B + C + BUFFER, BUFFER_SIZE, [
        DummyRecord(start=len(BUFFER + A + B), stop=len(BUFFER + A + B + C), op_type='INV', svid='dupINVdup'), # ABC -> ABc
        DummyRecord(start=len(BUFFER + A + B), stop=len(BUFFER + A + B + C), op_type='COPY-PASTE', svid='dupINVdup', target=len(BUFFER + A + B + C), insord=2), # ABC -> ABcC
        DummyRecord(start=len(BUFFER), stop=len(BUFFER + A), op_type='COPYinv-PASTE', svid='dupINVdup', target=len(BUFFER + A + B + C), insord=1), # ABC -> ABcaC
        DummyRecord(start=len(BUFFER + A), stop=len(BUFFER + A + B), op_type='CUTinv-PASTE', svid='dupINVdup', target=len(BUFFER + A + B + C), insord=0)], # ABC -> AcbaC
        [BUFFER + A + c + b + a + C + BUFFER]),
    GrammarCase('delINVdel', BUFFER + A + B + C + BUFFER, BUFFER_SIZE, [
        DummyRecord(start=len(BUFFER + A + B), stop=len(BUFFER + A + B + C), op_type='CUT', svid='delINVdel'), # ABC -> AB
        DummyRecord(start=len(BUFFER + A), stop=len(BUFFER + A + B), op_type='INV', svid='delINVdel'), # ABC -> Ab
        DummyRecord(start=len(BUFFER), stop=len(BUFFER + A), op_type='CUT', svid='delINVdel')], # ABC -> b
        [BUFFER + b + BUFFER]),
    GrammarCase('delINVdup', BUFFER + A + B + C + BUFFER, BUFFER_SIZE, [
        DummyRecord(start=len(BUFFER + A + B), stop=len(BUFFER + A + B + C), op_type='COPY-PASTE', svid='delINVdup', target=len(BUFFER + A + B + C), insord=1), # ABC -> ABCC
        DummyRecord(start=len(BUFFER + A), stop=len(BUFFER + A + B), op_type='CUTinv-PASTE', svid='delINVdup', target=len(BUFFER + A + B + C), insord=0), # ABC-> ABCbC ->ABcbC -> AcbC
        DummyRecord(start=len(BUFFER + A + B), stop=len(BUFFER + A + B + C), op_type='INV', svid='delINVdup'), 
        DummyRecord(start=len(BUFFER), stop=len(BUFFER + A), op_type='CUT', svid='delINVdup')], # ABC -> cbC
        [BUFFER + c + b + C + BUFFER]),
    GrammarCase('dupINVdel', BUFFER + A + B + C + BUFFER, BUFFER_SIZE, [
        DummyRecord(start=len(BUFFER + A + B), stop=len(BUFFER + A + B + C), op_type='CUT', svid='dupINVdel'), # ABC -> AB
        DummyRecord(start=len(BUFFER), stop=len(BUFFER + A), op_type='COPYinv-PASTE', svid='dupINVdel', target=len(BUFFER + A + B), insord=0), # ABC -> ABa
        DummyRecord(start=len(BUFFER + A), stop=len(BUFFER + A + B), op_type='INV', svid='dupINVdel')], # ABC -> Aba
        [BUFFER + A + b + a + BUFFER]),

    # dispersed, twice each: SHORT_DISP keeps one window, LONG_DISP splits into two
    GrammarCase('dDUP-short', BUFFER + A + SHORT_DISP + BUFFER, BUFFER_SIZE, [
        DummyRecord(start=len(BUFFER), stop=len(BUFFER + A), sv_type='dDUP', target=len(BUFFER + A + SHORT_DISP), insord=0)],
        [BUFFER + A + SHORT_DISP + A + BUFFER]),
    GrammarCase('dDUP-long', BUFFER + A + LONG_DISP + BUFFER, BUFFER_SIZE, [
        DummyRecord(start=len(BUFFER), stop=len(BUFFER + A), sv_type='dDUP', target=len(BUFFER + A + LONG_DISP), insord=0)],
        [BUFFER + A + DISP_HEAD,
         DISP_TAIL + A + BUFFER]),
    GrammarCase('INV_dDUP-short', BUFFER + A + SHORT_DISP + BUFFER, BUFFER_SIZE, [
        DummyRecord(start=len(BUFFER), stop=len(BUFFER + A), sv_type='INV_dDUP', target=len(BUFFER + A + SHORT_DISP), insord=0)],
        [BUFFER + A + SHORT_DISP + a + BUFFER]),
    GrammarCase('INV_dDUP-long', BUFFER + A + LONG_DISP + BUFFER, BUFFER_SIZE, [
        DummyRecord(start=len(BUFFER), stop=len(BUFFER + A), sv_type='INV_dDUP', target=len(BUFFER + A + LONG_DISP), insord=0)],
        [BUFFER + A + DISP_HEAD,
         DISP_TAIL + a + BUFFER]),
    GrammarCase('dDUP_INV-short', BUFFER + A + SHORT_DISP + BUFFER, BUFFER_SIZE, [
        DummyRecord(start=len(BUFFER), stop=len(BUFFER + A), op_type='INV', svid='dDUP_INV-short'),
        DummyRecord(start=len(BUFFER), stop=len(BUFFER + A), op_type='COPYinv-PASTE', svid='dDUP_INV-short', target=len(BUFFER + A + SHORT_DISP), insord=0)],
        [BUFFER + a + SHORT_DISP + a + BUFFER]),
    GrammarCase('dDUP_INV-long', BUFFER + A + LONG_DISP + BUFFER, BUFFER_SIZE, [
        DummyRecord(start=len(BUFFER), stop=len(BUFFER + A), op_type='INV', svid='dDUP_INV-long'),
        DummyRecord(start=len(BUFFER), stop=len(BUFFER + A), op_type='COPYinv-PASTE', svid='dDUP_INV-long', target=len(BUFFER + A + LONG_DISP), insord=0)],
        [BUFFER + a + DISP_HEAD,
         DISP_TAIL + a + BUFFER]),
    GrammarCase('INV_nrTRA-short', BUFFER + A + SHORT_DISP + BUFFER, BUFFER_SIZE, [
        DummyRecord(start=len(BUFFER), stop=len(BUFFER + A), sv_type='INV_nrTRA', target=len(BUFFER + A + SHORT_DISP), insord=0)],
        [BUFFER + SHORT_DISP + a + BUFFER]),
    GrammarCase('INV_nrTRA-long', BUFFER + A + LONG_DISP + BUFFER, BUFFER_SIZE, [
        DummyRecord(start=len(BUFFER), stop=len(BUFFER + A), sv_type='INV_nrTRA', target=len(BUFFER + A + LONG_DISP), insord=0)],
        [BUFFER + DISP_HEAD,
         DISP_TAIL + a + BUFFER]),
    GrammarCase('nrTRA-short', BUFFER + A + SHORT_DISP + BUFFER, BUFFER_SIZE, [
        DummyRecord(start=len(BUFFER), stop=len(BUFFER + A), sv_type='nrTRA', target=len(BUFFER + A + SHORT_DISP), insord=0)],
        [BUFFER + SHORT_DISP + A + BUFFER]),
    GrammarCase('nrTRA-long', BUFFER + A + LONG_DISP + BUFFER, BUFFER_SIZE, [
        DummyRecord(start=len(BUFFER), stop=len(BUFFER + A), sv_type='nrTRA', target=len(BUFFER + A + LONG_DISP), insord=0)],
        [BUFFER + DISP_HEAD,
         DISP_TAIL + A + BUFFER]),
    GrammarCase('rTRA-short', BUFFER + A + SHORT_DISP + B + BUFFER, BUFFER_SIZE, [
        DummyRecord(start=len(BUFFER), stop=len(BUFFER + A), op_type='CUT-PASTE', svid='rTRA-short', target=len(BUFFER + A + SHORT_DISP), insord=1), # 
        DummyRecord(start=len(BUFFER + A + SHORT_DISP), stop=len(BUFFER + A + SHORT_DISP + B), op_type='CUT-PASTE', svid='rTRA-short', target=len(BUFFER), insord=0)],
        [BUFFER + B + SHORT_DISP + A + BUFFER]),
    GrammarCase('rTRA-long', BUFFER + A + LONG_DISP + B + BUFFER, BUFFER_SIZE, [
        DummyRecord(start=len(BUFFER), stop=len(BUFFER + A), op_type='CUT-PASTE', svid='rTRA-long', target=len(BUFFER + A + LONG_DISP), insord=1),
        DummyRecord(start=len(BUFFER + A + LONG_DISP), stop=len(BUFFER + A + LONG_DISP + B), op_type='CUT-PASTE', svid='rTRA-long', target=len(BUFFER), insord=0)],
        [BUFFER + B + DISP_HEAD,
         DISP_TAIL + A + BUFFER]),
    GrammarCase('INV_rTRA-short', BUFFER + A + SHORT_DISP + B + BUFFER, BUFFER_SIZE, [
        DummyRecord(start=len(BUFFER), stop=len(BUFFER + A), op_type='CUTinv-PASTE', svid='INV_rTRA-short', target=len(BUFFER + A + SHORT_DISP), insord=1),
        DummyRecord(start=len(BUFFER + A + SHORT_DISP), stop=len(BUFFER + A + SHORT_DISP + B), op_type='CUTinv-PASTE', svid='INV_rTRA-short', target=len(BUFFER), insord=0)],
        [BUFFER + b + SHORT_DISP + a + BUFFER]),
    GrammarCase('INV_rTRA-long', BUFFER + A + LONG_DISP + B + BUFFER, BUFFER_SIZE, [
        DummyRecord(start=len(BUFFER), stop=len(BUFFER + A), op_type='CUTinv-PASTE', svid='INV_rTRA-long', target=len(BUFFER + A + LONG_DISP), insord=1),
        DummyRecord(start=len(BUFFER + A + LONG_DISP), stop=len(BUFFER + A + LONG_DISP + B), op_type='CUTinv-PASTE', svid='INV_rTRA-long', target=len(BUFFER), insord=0)],
        [BUFFER + b + DISP_HEAD,
         DISP_TAIL + a + BUFFER]),
    GrammarCase('INS_iDEL-short', BUFFER + A + SHORT_DISP + B + BUFFER, BUFFER_SIZE, [
        DummyRecord(start=len(BUFFER + A + SHORT_DISP), stop=len(BUFFER + A + SHORT_DISP + B), op_type='CUT', svid='INS_iDEL-short'),
        DummyRecord(start=len(BUFFER), stop=len(BUFFER + A), op_type='CUT-PASTE', svid='INS_iDEL-short', target=len(BUFFER + A + SHORT_DISP), insord=0)],
        [BUFFER + SHORT_DISP + A + BUFFER]),
    GrammarCase('INS_iDEL-long', BUFFER + A + LONG_DISP + B + BUFFER, BUFFER_SIZE, [
        DummyRecord(start=len(BUFFER), stop=len(BUFFER + A), op_type='CUT-PASTE', svid='INS_iDEL-long', target=len(BUFFER + A + LONG_DISP), insord=0),
        DummyRecord(start=len(BUFFER + A + LONG_DISP), stop=len(BUFFER + A + LONG_DISP + B), op_type='CUT', svid='INS_iDEL-long')],
        [BUFFER + DISP_HEAD,
         DISP_TAIL + A + BUFFER]),
    GrammarCase('dDUP_iDEL-short', BUFFER + A + SHORT_DISP + B + BUFFER, BUFFER_SIZE, [
        DummyRecord(start=len(BUFFER), stop=len(BUFFER + A), op_type='COPY-PASTE', svid='dDUP_iDEL-short', target=len(BUFFER + A + SHORT_DISP), insord=0),
        DummyRecord(start=len(BUFFER + A + SHORT_DISP), stop=len(BUFFER + A + SHORT_DISP + B), op_type='CUT', svid='dDUP_iDEL-short')],
        [BUFFER + A + SHORT_DISP + A + BUFFER]),
    GrammarCase('dDUP_iDEL-long', BUFFER + A + LONG_DISP + B + BUFFER, BUFFER_SIZE, [
        DummyRecord(start=len(BUFFER + A + LONG_DISP), stop=len(BUFFER + A + LONG_DISP + B), op_type='CUT', svid='dDUP_iDEL-long'),
        DummyRecord(start=len(BUFFER), stop=len(BUFFER + A), op_type='COPY-PASTE', svid='dDUP_iDEL-long', target=len(BUFFER + A + LONG_DISP), insord=0)],
        [BUFFER + A + DISP_HEAD,
         DISP_TAIL + A + BUFFER]),
]


@pytest.mark.parametrize('reference,buffer_size,records,expected_sequences',
                         [case[1:] for case in GRAMMAR_CASES],
                         ids=[case.name for case in GRAMMAR_CASES])
def test_reconstruction_matches_grammar(reference, buffer_size, records,
                                        expected_sequences):
    """Records built directly from each type's known grammar, grouped as the callset loader does."""
    reference_bytes = {SYNTHETIC_CHROM: reference.encode()}

    grouped_records = group_records_by_id(records)
    # every case is one SV, whether by shared SVID or as a simple variant
    assert len(grouped_records) == 1
    (_, sv_records), = grouped_records.items()

    queries = construct_queries(sv_records, buffer_size, reference_bytes)

    assert [query.sequence for query in queries] == expected_sequences
