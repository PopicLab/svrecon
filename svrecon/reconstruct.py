"""Alt-allele reconstruction: transform the reference window into the SV's alt allele and record,
for each piece of the result, a ``(start, end)`` segment (in alt-allele coordinates) for scoring's
``validate_segments_from_cigar`` to check locally.

Segments cover the whole allele: inserted and inverted pieces, and the untouched reference runs
between them (including the left/right context buffer). A deletion simply removes its segment,
joining its neighbors directly. Segments are independent, non-overlapping ranges; an overlap
(e.g. from imprecise caller coordinates) raises ValueError rather than being silently resolved.
"""
import bisect
import logging
from collections import namedtuple
from dataclasses import dataclass
from typing import Dict, List, Tuple

from pysam import VariantRecord

from svrecon.util import get_start_stop, reverse_complement

logger = logging.getLogger(__name__)


@dataclass
class QueryReconSubsequence:
    chrom: str
    svtype: str
    svid: str
    sequence: str
    length: int
    ref_start: int
    ref_end: int
    segments: List[Tuple[int, int]] # indicator of segments of interest, 0 indexed to the start of the subsequence,
    ref_segments: List[Tuple[int, int]] # same pieces, in absolute reference coordinates (parallel to segments)

    def __len__(self):
        return self.length


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
    invert: bool = False
    insord: int = -1

    def modify_segments(self, segments: List["_Segment"]) -> None:
        length = self.ref_end - self.ref_start
        idx = bisect.bisect_right(segments, self.op_start, key=lambda seg: seg.alt_end)
        segments.insert(idx, _Segment(ref_start=self.ref_start, ref_end=self.ref_end, alt_start=self.op_start,
                                      alt_end=self.op_start + length, invert=self.invert))
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

def get_operations(records: list[VariantRecord]) -> list[_Operation]:
    operations: list[_Operation] = []

    for record in records:
        start, stop = get_start_stop(record)
        target = record.info.get('TARGET', stop) # NOTE: is there a +1 here?
        insord = record.info.get('INSORD', -1)

        if record.info['OP_TYPE'] == 'CUT' or record.info['SVTYPE'] == 'DEL':
            operations.append(_Delete(start, start, stop))
        elif record.info['OP_TYPE'] == 'INV' or record.info['SVTYPE'] == 'INV':
            operations.append(_Invert(start, start, stop))
        elif record.info['OP_TYPE'] == 'DUP' or record.info['SVTYPE'] == 'DUP':
            operations.append(_Insert(stop, start, stop))
        elif record.info['OP_TYPE'] == 'COPY-PASTE' or record.info['SVTYPE'] == 'dDUP':
            operations.append(_Insert(target, start, stop, insord=insord))
        elif record.info['OP_TYPE'] == 'CUT-PASTE' or record.info['SVTYPE'] == 'nrTRA':
            operations.append(_Delete(start, start, stop))
            operations.append(_Insert(target, start, stop, insord=insord))
        elif record.info['OP_TYPE'] == 'COPYinv-PASTE' or record.info['SVTYPE'] == 'INV_dDUP':
            operations.append(_Insert(target, start, stop, invert=True, insord=insord))
        elif record.info['OP_TYPE'] == 'CUTinv-PASTE' or record.info['SVTYPE'] == 'INV_nrTRA':
            operations.append(_Delete(start, start, stop))
            operations.append(_Insert(target, start, stop, invert=True, insord=insord))
        else:
            logger.warning(f'Unknown OP_TYPE: {record.info["OP_TYPE"]}')

    return operations

def create_starting_segments(records: list[VariantRecord], buffer: int) -> list[_Segment]:
    """
    Creates buffer segments around boundaries of contiguous records
    """
    starts_stops = sorted({get_start_stop(rec) for rec in records})
    min_start = starts_stops[0][0]
    max_stop = max(stop for _, stop in starts_stops)

    segments = [
        _Segment(max(0, min_start - buffer), min_start, max(0, min_start - buffer), min_start, False),
        _Segment(max_stop, max_stop + buffer, max_stop, max_stop + buffer, False),
    ]
    segments.extend(_Segment(start, stop, start, stop, False) for start, stop in starts_stops)

    targets = {rec.info['TARGET'] for rec in records if 'TARGET' in rec.info}
    for target in sorted(targets):
        if target < min_start or target > max_stop:
            segments.append(_Segment(max(0, target - buffer), target, max(0, target - buffer), target, False))
            segments.append(_Segment(target, target + buffer, target, target + buffer, False))

    segments.sort(key=lambda seg: seg.alt_start)
    for prev, seg in zip(segments, segments[1:]):
        if seg.alt_start < prev.alt_end:
            raise ValueError(f'starting segment [{seg.alt_start},{seg.alt_end}) overlaps preceding segment '
                             f'[{prev.alt_start},{prev.alt_end}) by {prev.alt_end - seg.alt_start} bp')

    return segments

def simulate_subsequences(records: List[VariantRecord], buffer: int, ref: Dict[str, bytearray]) -> List[QueryReconSubsequence]:
    """
    Given records of operations, recreate the resulting subsequences and metadata 
    """

    # Descending by op_start type order (Invert, Delete, Insert): INSORD
    _OP_TYPE_ORDER = {_Invert: 0, _Delete: 1, _Insert: 2}
    operations: list[_Operation] = sorted(
        get_operations(records),
        key=lambda op: (-op.op_start, _OP_TYPE_ORDER[type(op)],
                         -op.insord if isinstance(op, _Insert) else 0))
    segments: list[_Segment] = create_starting_segments(records, buffer)

    for operation in operations:
        operation.modify_segments(segments)

    sv_type = records[0].info['SVTYPE']
    svid = records[0].info['SVID']
    chrom = records[0].chrom

    # Split into maximal contiguous runs -- a gap marks a separate, independent window
    # (e.g. a dispersed duplication's distant destination) -- and materialize each run's
    # bases from ref into its own QueryReconSubsequence.
    runs: List[List[_Segment]] = []
    for seg in segments:
        if runs and seg.alt_start == runs[-1][-1].alt_end:
            runs[-1].append(seg)
        else:
            runs.append([seg])

    subsequences = []
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
        subsequences.append(QueryReconSubsequence(
            chrom=chrom, svtype=sv_type, svid=svid, sequence=sequence, length=len(sequence),
            ref_start=run[0].ref_start, ref_end=run[-1].ref_end, segments=covering_segments,
            ref_segments=covering_ref_segments))

    return subsequences