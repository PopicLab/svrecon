"""Sequence reconstruction"""
import bisect
from dataclasses import dataclass
from typing import Dict, List, Tuple
import string

from pysam import VariantRecord

from svrecon.utils import get_start_stop, reverse_complement


@dataclass
class Query:
    chrom: str
    svtype: str
    svid: str
    grammar: str
    sequence: str
    ref_start: int
    ref_end: int
    recon_segments: List[Tuple[int, int]] # indicator of segments of interest, 0 indexed to the start of the subsequence,
    ref_segments: List[Tuple[int, int]] # same pieces, in absolute reference coordinates (parallel to segments)
    ref_sequence: str # unmodified reference over [ref_start, ref_end), for plotting and diagnostics
    buffer: int # bp length of reference context included on each side of SV regions

    def __len__(self):
        return len(self.sequence)

@dataclass
class _Segment:
    """
    Contiguous [start, end) intervals.
    ref* - positions of the original sequence
    alt* - positions of the new sequence, same coordinate system as ref
    invert - inverted sequence flag
    """
    ref_start: int
    ref_end: int
    alt_start: int
    alt_end: int
    invert: bool
    symbol: str

@dataclass
class _Operation:
    """
    Performs operation on a list of segments at position op_start,
    """
    op_start: int
    ref_start: int
    ref_end: int

    def modify_segments(self, segments: List[_Segment]) -> None:
        """
        Modifies segments in place, maintaining [alt_start, alt_end) sort, which represents a sequence.
        Upon modification, adjust all following segments alt_start and alt_end by the length of the change.
        """
        raise NotImplementedError("Implement in subclass")

@dataclass
class _Insert(_Operation):
    symbol: str
    invert: bool = False
    insord: int = -1

    def modify_segments(self, segments: List["_Segment"]) -> None:
        length = self.ref_end - self.ref_start
        idx = bisect.bisect_right(segments, self.op_start, key=lambda seg: seg.alt_end)
        symbol = self.symbol.lower() if self.invert else self.symbol
        segments.insert(idx, _Segment(ref_start=self.ref_start, ref_end=self.ref_end, alt_start=self.op_start,
                                      alt_end=self.op_start + length, invert=self.invert, symbol=symbol))
        for seg in segments[idx + 1:]:
            seg.alt_start += length
            seg.alt_end += length

@dataclass
class _Delete(_Operation):
    def modify_segments(self, segments: List["_Segment"]) -> None:
        length = self.ref_end - self.ref_start
        idx = bisect.bisect_right(segments, self.op_start, key=lambda seg: seg.alt_end)
        seg = segments[idx]
        if seg.alt_start != self.op_start or seg.alt_end != self.op_start + length:
            raise ValueError(f'Delete at alt position {self.op_start} (length {length}) does not land exactly '
                             f'on an existing segment boundary (found [{seg.alt_start},{seg.alt_end}))')
        segments.pop(idx)
        for later in segments[idx:]:
            later.alt_start -= length
            later.alt_end -= length

@dataclass
class _Invert(_Operation):
    def modify_segments(self, segments: List["_Segment"]) -> None:
        length = self.ref_end - self.ref_start
        idx = bisect.bisect_right(segments, self.op_start, key=lambda seg: seg.alt_end)
        seg = segments[idx]
        if seg.alt_start != self.op_start or seg.alt_end != self.op_start + length:
            raise ValueError(f'Invert at alt position {self.op_start} (length {length}) does not land exactly '
                             f'on an existing segment boundary (found [{seg.alt_start},{seg.alt_end}))')
        seg.invert = not seg.invert
        seg.symbol = seg.symbol.swapcase()

def get_operations_from_records(records: list[VariantRecord]) -> list[_Operation]:
    operations: list[_Operation] = []
    starts_stops = sorted(set(get_start_stop(rec) for rec in records))
    start_stops_to_sym: dict[Tuple, str] = {interval: sym for interval, sym in zip(starts_stops, string.ascii_uppercase)}

    for record in records:
        start, stop = interval = get_start_stop(record)
        target = record.info.get('TARGET', stop) # already 0-based: TARGET == start + SVLEN for tandem pastes
        insord = record.info.get('INSORD', -1)
        # if single record SV, read off SVTYPE. otherwise, read OP_TYPE
        op = record.info['SVTYPE'] if len(records) == 1 else record.info['OP_TYPE']
        symbol = start_stops_to_sym[interval]

        if op in ('CUT', 'DEL'):
            operations.append(_Delete(start, start, stop))
        elif op == 'INV':
            operations.append(_Invert(start, start, stop))
        elif op == 'DUP':
            operations.append(_Insert(stop, start, stop, symbol))
        elif op in ('COPY-PASTE', 'dDUP'):
            operations.append(_Insert(target, start, stop, symbol=symbol, insord=insord))
        elif op in ('CUT-PASTE', 'nrTRA'):
            operations.append(_Delete(start, start, stop))
            operations.append(_Insert(target, start, stop, insord=insord, symbol=symbol))
        elif op in ('COPYinv-PASTE', 'INV_dDUP', 'INV_DUP'):
            operations.append(_Insert(target, start, stop, invert=True, insord=insord, symbol=symbol))
        elif op in ('CUTinv-PASTE', 'INV_nrTRA'):
            operations.append(_Delete(start, start, stop))
            operations.append(_Insert(target, start, stop, invert=True, insord=insord, symbol=symbol))
        else:
            raise ValueError(f'Unknown operation type: {op}')

    return operations

def create_starting_segments(records: list[VariantRecord], buffer: int, ref: Dict[str, bytearray]) -> list[_Segment]:
    """
    Creates buffer segments around boundaries of contiguous records, clamped to [0, chromosome length)
    """
    starts_stops = sorted(set(get_start_stop(rec) for rec in records))
    start_stops_to_sym: dict[Tuple, str] = {interval: sym for interval, sym in zip(starts_stops, string.ascii_uppercase)}

    chrom_len = len(ref[records[0].chrom])
    targets = sorted({rec.info['TARGET'] for rec in records if 'TARGET' in rec.info})
    breakpoints = sorted({pos for start, stop in starts_stops for pos in (start, stop)} | set(targets))
    clamp = lambda pos: min(max(pos, 0), chrom_len)

    # Preliminary buffers surrounding source span and targets
    windows = [(clamp(start - buffer), clamp(stop + buffer)) for start, stop in starts_stops]
    windows += [(clamp(target - buffer), clamp(target + buffer)) for target in targets]

    # Merge overlapping sources & targets with buffers
    merged: List[Tuple[int, int]] = []
    for window_start, window_stop in sorted(windows):
        if merged and window_start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], window_stop))
        else:
            merged.append((window_start, window_stop))

    # Cut each merged buffer at the breakpoints inside it
    segments: List[_Segment] = []
    for window_start, window_stop in merged:
        cuts = [window_start]
        for breakpoint_pos in breakpoints:
            if window_start < breakpoint_pos < window_stop:
                cuts.append(breakpoint_pos)
        cuts.append(window_stop)

        for segment_start, segment_stop in zip(cuts, cuts[1:]):
            segments.append(_Segment(ref_start=segment_start, ref_end=segment_stop,
                                     alt_start=segment_start, alt_end=segment_stop, invert=False, 
                                     symbol=start_stops_to_sym.get((segment_start, segment_stop), '~'))) # assign ~ as buffer regions

    return segments

def construct_queries(svid: str, records: List[VariantRecord], buffer: int,
                      ref: Dict[str, bytearray]) -> List[Query]:
    """
    Given records of operations, recreate the resulting subsequences and metadata 
    """

    # Descending by op_start type order (Invert, Delete, Insert): INSORD
    _OP_TYPE_ORDER = {_Invert: 0, _Delete: 1, _Insert: 2}
    operations: list[_Operation] = sorted(
        get_operations_from_records(records),
        key=lambda op: (-op.op_start, _OP_TYPE_ORDER[type(op)],
                         -op.insord if isinstance(op, _Insert) else 0))
    segments: list[_Segment] = create_starting_segments(records, buffer, ref)
    initial_grammar = ''.join(seg.symbol for seg in segments) # before operations mutate segments; same for every query of this SV

    for operation in operations:
        operation.modify_segments(segments)

    sv_type = records[0].info['SVTYPE']
    chrom = records[0].chrom

    # Split into maximal contiguous runs -- a gap marks a separate, independent window
    runs: List[List[_Segment]] = []
    for seg in segments:
        if runs and seg.alt_start == runs[-1][-1].alt_end:
            runs[-1].append(seg)
        else:
            runs.append([seg])

    queries = []
    for run in runs:
        offset = run[0].alt_start
        pieces = []
        covering_segments = []
        covering_ref_segments = []
        for seg in run:
            clip = ref[chrom][seg.ref_start:seg.ref_end].decode('ascii')
            if seg.invert:
                clip = reverse_complement(clip)
            pieces.append(clip)
            covering_segments.append((seg.alt_start - offset, seg.alt_end - offset))
            covering_ref_segments.append((seg.ref_start, seg.ref_end))

        sequence = ''.join(pieces)
        resultant_grammar = ''.join(seg.symbol for seg in run) # include buffers for debugging
        grammar = f'{initial_grammar}->{resultant_grammar}'
        query_ref_start, query_ref_end = run[0].ref_start, run[-1].ref_end
        queries.append(Query(
            chrom=chrom, svtype=sv_type, svid=svid, grammar=grammar, sequence=sequence,
            ref_start=query_ref_start, ref_end=query_ref_end, recon_segments=covering_segments,
            ref_segments=covering_ref_segments,
            ref_sequence=ref[chrom][query_ref_start:query_ref_end].decode('ascii'),
            buffer=buffer))

    return queries