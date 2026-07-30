"""Alt-allele reconstruction: transform the reference window into the SV's alt allele and record,
for each piece of the result, a ``(start, end)`` segment (in alt-allele coordinates) for scoring's
``validate_segments_from_cigar`` to check locally.

Segments cover the whole allele: inserted and inverted pieces, and the untouched reference runs
between them (including the left/right context buffer). A deletion is simply a zero-length segment
at its join. Segments are independent ranges, so pieces that overlap because of imprecise caller
coordinates are representable and reconciled by clipping the earlier piece (a clip larger than
``SEGMENT_CLIP_WARN_THRESHOLD`` bp is logged, since that signals a real coordinate error).
"""
import logging
from dataclasses import dataclass
from typing import Dict, List, Tuple

from pysam import VariantRecord

from svrecon.util import get_start_stop, reverse_complement

logger = logging.getLogger(__name__)

# A reconstructed piece clipped by more than this many bp to resolve an overlap probably reflects a
# real coordinate error rather than 1-2 bp of caller rounding -> warn.
SEGMENT_CLIP_WARN_THRESHOLD = 10


@dataclass
class QueryReconSubsequence:
    chrom: str
    svtype: str
    svid: str
    sequence: str
    location: int
    length: int
    ref_start: int
    ref_end: int
    segments: List[Tuple[int, int]]

    def __len__(self):
        return self.length


@dataclass
class _Segment:
    """A piece of the allele while it is being built; coordinates are mutable so operations can
    shift them in-place and overlaps can be clipped. A deletion is just a segment whose bases become ''
    placeholders, so it collapses to zero length once the placeholders are dropped."""
    start: int
    end: int


def simulate_subsequences(records: List[VariantRecord], buffer: int, ref: Dict[str, bytearray]) -> List[QueryReconSubsequence]:

    sv_type = records[0].info['SVTYPE']
    svid = records[0].info['SVID']
    chrom = records[0].chrom

    target_records = []
    in_place_records = []

    for rec in records:
        if rec.info['SVTYPE'] == 'DUP':
            rec.info['TARGET'] = rec.stop

        if 'TARGET' in rec.info:
            target_records.append(rec)
        else:
            in_place_records.append(rec)

    if len(target_records) > 1:
        target_records = sorted(target_records, key=lambda rec: (rec.info['TARGET'], rec.info.get('INSORD', 0)),
                                reverse=True)

    records = in_place_records + target_records

    offset = max(0, min([rec.start for rec in records]) - buffer)
    sequence_end = max([rec.stop for rec in records]) + buffer

    target_offset = max(0, min([rec.info['TARGET'] for rec in target_records]) - buffer) if target_records else 0
    target_sequence_end = max([rec.info['TARGET'] for rec in target_records]) + buffer if target_records else 0

    merged_target_sequence = offset <= target_offset + buffer and target_sequence_end - buffer <= sequence_end
    if merged_target_sequence:
        target_offset = offset

    orig_sequence = list(ref[chrom][offset:sequence_end].decode('ascii'))
    new_sequence = orig_sequence.copy()
    # Per-window segments, in that window's LIST-index coordinates. A deletion keeps its bases as ''
    # placeholders (so later ops' reference coordinates stay valid) and collapses to a zero-length
    # segment once the placeholders are dropped in _finalize_window.
    main_segments: List[_Segment] = []

    if not merged_target_sequence:
        new_target_sequence = list(ref[chrom][target_offset:target_sequence_end].decode('ascii'))
        target_segments: List[_Segment] = []

    delete_placeholder = ''
    queries = []

    def shift_after(segments, position, amount):
        """An insertion of `amount` bases at `position` pushes everything at/after it to the right."""
        for segment in segments:
            if segment.start >= position:
                segment.start += amount
            if segment.end > position:
                segment.end += amount

    # Pass 1: every operation's SOURCE region is a segment, in original (pre-transform) coordinates.
    # Creating them up front means the pastes in pass 2 shift them along with everything else after
    # the insertion point, so a source that ends up after an insertion still lands at the right place.
    for rec in records:
        if rec.info.get('TARGET_CHROM', chrom) != chrom:
            logger.warning(
                f'Skipping record: {rec.id}-{sv_type} because interchromosome target. Interchromosome checks not implemented yet.')
            return []
        start, stop = get_start_stop(rec)
        main_segments.append(_Segment(start - offset, stop - offset))

    # Pass 2: apply the transformations. In-place ops (DEL/INV) only rewrite bases -- their segment
    # already exists. A paste inserts its copy, shifts every segment after the insertion (sources
    # included), and records the copy.
    for rec in records:
        start, stop = get_start_stop(rec)
        target = rec.info.get('TARGET', stop + 1) - target_offset
        start -= offset
        stop -= offset

        if rec.info['OP_TYPE'] == 'CUT' or rec.info['SVTYPE'] == 'DEL':
            new_sequence[start:stop] = [delete_placeholder] * (stop - start)
        elif rec.info['OP_TYPE'] == 'INV' or rec.info['SVTYPE'] == 'INV':
            new_sequence[start:stop] = reverse_complement(orig_sequence[start:stop])
        elif rec.info['OP_TYPE'] == 'COPY-PASTE' or rec.info['SVTYPE'] in ['DUP', 'dDUP']:
            clip = orig_sequence[start:stop]
            if merged_target_sequence:
                new_sequence = new_sequence[:target] + clip + new_sequence[target:]
                shift_after(main_segments, target, len(clip))
                main_segments.append(_Segment(target, target + len(clip)))
            else:
                new_target_sequence = new_target_sequence[:target] + clip + new_target_sequence[target:]
                shift_after(target_segments, target, len(clip))
                target_segments.append(_Segment(target, target + len(clip)))
        elif rec.info['OP_TYPE'] in ['CUT-PASTE'] or rec.info['SVTYPE'] == 'nrTRA':
            clip = orig_sequence[start:stop]
            new_sequence[start:stop] = [delete_placeholder] * (stop - start)
            if merged_target_sequence:
                new_sequence = new_sequence[:target] + clip + new_sequence[target:]
                shift_after(main_segments, target, len(clip))
                main_segments.append(_Segment(target, target + len(clip)))
            else:
                new_target_sequence = new_target_sequence[:target] + clip + new_target_sequence[target:]
                shift_after(target_segments, target, len(clip))
                target_segments.append(_Segment(target, target + len(clip)))
        elif rec.info['OP_TYPE'] == 'COPYinv-PASTE' or rec.info['SVTYPE'] in ['INV_dDUP']:
            clip = list(reverse_complement(orig_sequence[start:stop]))
            if merged_target_sequence:
                new_sequence = new_sequence[:target] + clip + new_sequence[target:]
                shift_after(main_segments, target, len(clip))
                main_segments.append(_Segment(target, target + len(clip)))
            else:
                new_target_sequence = new_target_sequence[:target] + clip + new_target_sequence[target:]
                shift_after(target_segments, target, len(clip))
                target_segments.append(_Segment(target, target + len(clip)))
        elif rec.info['OP_TYPE'] == 'CUTinv-PASTE' or rec.info['SVTYPE'] in ['INV_nrTRA']:
            clip = list(reverse_complement(orig_sequence[start:stop]))
            new_sequence[start:stop] = [delete_placeholder] * (stop - start)
            if merged_target_sequence:
                new_sequence = new_sequence[:target] + clip + new_sequence[target:]
                shift_after(main_segments, target, len(clip))
                main_segments.append(_Segment(target, target + len(clip)))
            else:
                new_target_sequence = new_target_sequence[:target] + clip + new_target_sequence[target:]
                shift_after(target_segments, target, len(clip))
                target_segments.append(_Segment(target, target + len(clip)))
        else:
            logger.warning(f'Unknown OP_TYPE: {rec.info["OP_TYPE"]}')

    queries.append(_finalize_window(new_sequence, main_segments, offset, sequence_end, svid, chrom, sv_type))
    if not merged_target_sequence:
        queries.append(_finalize_window(new_target_sequence, target_segments, target_offset, target_sequence_end,
                                        svid, chrom, sv_type))

    return queries


def _finalize_window(window_chars, segments, offset, ref_end, svid, chrom, sv_type) -> QueryReconSubsequence:
    """Assemble one transformed window into a QueryReconSubsequence: the reconstructed pieces plus
    the SV's two outer flanks, in allele-string coordinates."""
    # Join the transformed characters into the allele and translate each segment from LIST indices to
    # allele-string positions. Deletion placeholders ('') contribute no characters, so a deleted run
    # collapses to a zero-length segment at its join.
    string_position_at = []          # string_position_at[list_index] -> allele-string offset
    string_position = 0
    for base in window_chars:
        string_position_at.append(string_position)
        string_position += len(base)  # '' placeholder -> 0
    string_position_at.append(string_position)  # end sentinel
    sequence = ''.join(window_chars)
    segments = [_Segment(string_position_at[s.start], string_position_at[s.end]) for s in segments]

    # Reconcile overlaps: pieces should be disjoint, but imprecise caller coordinates can make two
    # overlap; keep the later-created piece's boundary and clip the earlier one.
    placed = []
    for segment in segments:  # in creation order
        for earlier in placed:
            _clip_out_overlap(earlier, segment, svid, sv_type)
        placed.append(segment)

    # Cover the allele with the pieces plus the SV's outer flanks. 
    pieces = [(seg.start, seg.end) for seg in placed]
    sv_start = min(start for start, _ in pieces)
    sv_end = max(end for _, end in pieces)
    flanks = [flank for flank in [(0, sv_start), (sv_end, len(sequence))] if flank[1] > flank[0]]
    covering_segments = sorted(set(pieces + flanks))

    return QueryReconSubsequence(
        chrom=chrom, svtype=sv_type, svid=svid, sequence=sequence,
        location=offset, length=len(sequence), ref_start=offset, ref_end=ref_end,
        segments=covering_segments)


def _clip_out_overlap(earlier: _Segment, newer: _Segment, svid, sv_type):
    """Trim `earlier` so it no longer overlaps `newer` (which owns the boundary)."""
    overlap = min(earlier.end, newer.end) - max(earlier.start, newer.start)
    if overlap <= 0:
        return
    if overlap > SEGMENT_CLIP_WARN_THRESHOLD:
        logger.warning(f'{svid}-{sv_type}: reconstructed segment [{earlier.start},{earlier.end}] clipped '
                       f'{overlap} bp by [{newer.start},{newer.end}] (imprecise caller coordinates?)')
    if earlier.start < newer.start:
        earlier.end = min(earlier.end, newer.start)
    else:
        earlier.start = max(earlier.start, newer.end)
