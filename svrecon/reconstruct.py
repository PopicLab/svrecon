"""Alt-allele reconstruction: build the changed subsequences (with junction masks) that scoring aligns against the assembly/reads."""
from dataclasses import dataclass
from itertools import groupby
from typing import Dict, List

from pysam import VariantRecord

from svrecon.util import get_start_stop, reverse_complement


@dataclass
class QueryReconSubsequence:
    chrom: str
    svtype: str
    svid: str
    sequence: str
    location: int
    length: int
    junctions: List[int]
    result_len: int


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
    changed_mask = [False] * len(new_sequence)
    junction_mask = [False] * len(new_sequence)

    if not merged_target_sequence:
        new_target_sequence = list(ref[chrom][target_offset:target_sequence_end].decode('ascii'))
        target_changed_mask = [False] * len(new_target_sequence)
        target_junction_mask = [False] * len(new_target_sequence)

    delete_placeholder = ''
    queries = []

    def mark_junction(mask, idx):
        if 0 <= idx < len(mask):
            mask[idx] = True

    for rec in records:
        start, stop = get_start_stop(rec)
        target = rec.info.get('TARGET', stop + 1) - target_offset
        start -= offset
        stop -= offset

        if rec.info.get('TARGET_CHROM', chrom) != chrom:
            logger.warning(
                f'Skipping record: {rec.id}-{sv_type} because interchromosome target. Interchromosome checks not implemented yet.')
            return []

        if rec.info['OP_TYPE'] == 'CUT' or rec.info['SVTYPE'] == 'DEL':
            new_sequence[start:stop] = [delete_placeholder] * (stop - start)
            changed_mask[start:stop] = [True] * (stop - start)
            mark_junction(junction_mask, start)
            mark_junction(junction_mask, stop - 1)
        elif rec.info['OP_TYPE'] == 'INV' or rec.info['SVTYPE'] == 'INV':
            new_sequence[start:stop] = reverse_complement(orig_sequence[start:stop])
            changed_mask[start:stop] = [True] * (stop - start)
            mark_junction(junction_mask, start)
            mark_junction(junction_mask, stop - 1)
        elif rec.info['OP_TYPE'] == 'COPY-PASTE' or rec.info['SVTYPE'] in ['DUP', 'dDUP']:
            clip = orig_sequence[start:stop]
            if merged_target_sequence:
                new_sequence = new_sequence[:target] + clip + new_sequence[target:]
                changed_mask = changed_mask[:target] + [True] * len(clip) + changed_mask[target:]
                junction_mask = junction_mask[:target] + [False] * len(clip) + junction_mask[target:]
                mark_junction(junction_mask, target)
                mark_junction(junction_mask, target + len(clip) - 1)
            else:
                new_target_sequence = new_target_sequence[:target] + clip + new_target_sequence[target:]
                target_changed_mask = target_changed_mask[:target] + [True] * len(clip) + target_changed_mask[target:]
                target_junction_mask = target_junction_mask[:target] + [False] * len(clip) + target_junction_mask[
                    target:]
                mark_junction(target_junction_mask, target)
                mark_junction(target_junction_mask, target + len(clip) - 1)
        elif rec.info['OP_TYPE'] in ['CUT-PASTE'] or rec.info['SVTYPE'] == 'nrTRA':
            clip = orig_sequence[start:stop]
            new_sequence[start:stop] = [delete_placeholder] * (stop - start)
            changed_mask[start:stop] = [True] * len(clip)
            mark_junction(junction_mask, start)
            mark_junction(junction_mask, stop - 1)

            if merged_target_sequence:
                new_sequence = new_sequence[:target] + clip + new_sequence[target:]
                changed_mask = changed_mask[:target] + [True] * len(clip) + changed_mask[target:]
                junction_mask = junction_mask[:target] + [False] * len(clip) + junction_mask[target:]
                mark_junction(junction_mask, target)
                mark_junction(junction_mask, target + len(clip) - 1)
            else:
                new_target_sequence = new_target_sequence[:target] + clip + new_target_sequence[target:]
                target_changed_mask = target_changed_mask[:target] + [True] * len(clip) + target_changed_mask[target:]
                target_junction_mask = target_junction_mask[:target] + [False] * len(clip) + target_junction_mask[
                    target:]
                mark_junction(target_junction_mask, target)
                mark_junction(target_junction_mask, target + len(clip) - 1)
        elif rec.info['OP_TYPE'] == 'COPYinv-PASTE' or rec.info['SVTYPE'] in ['INV_dDUP']:
            clip = list(reverse_complement(orig_sequence[start:stop]))

            if merged_target_sequence:
                new_sequence = new_sequence[:target] + clip + new_sequence[target:]
                changed_mask = changed_mask[:target] + [True] * len(clip) + changed_mask[target:]
                junction_mask = junction_mask[:target] + [False] * len(clip) + junction_mask[target:]
                mark_junction(junction_mask, target)
                mark_junction(junction_mask, target + len(clip) - 1)
            else:
                new_target_sequence = new_target_sequence[:target] + clip + new_target_sequence[target:]
                target_changed_mask = target_changed_mask[:target] + [True] * len(clip) + target_changed_mask[target:]
                target_junction_mask = target_junction_mask[:target] + [False] * len(clip) + target_junction_mask[
                    target:]
                mark_junction(target_junction_mask, target)
                mark_junction(target_junction_mask, target + len(clip) - 1)
        elif rec.info['OP_TYPE'] == 'CUTinv-PASTE' or rec.info['SVTYPE'] in ['INV_nrTRA']:
            clip = list(reverse_complement(orig_sequence[start:stop]))
            new_sequence[start:stop] = [delete_placeholder] * (stop - start)
            changed_mask[start:stop] = [True] * len(clip)
            mark_junction(junction_mask, start)
            mark_junction(junction_mask, stop - 1)

            if merged_target_sequence:
                new_sequence = new_sequence[:target] + clip + new_sequence[target:]
                changed_mask = changed_mask[:target] + [True] * len(clip) + changed_mask[target:]
                junction_mask = junction_mask[:target] + [False] * len(clip) + junction_mask[target:]
                mark_junction(junction_mask, target)
                mark_junction(junction_mask, target + len(clip) - 1)
            else:
                new_target_sequence = new_target_sequence[:target] + clip + new_target_sequence[target:]
                target_changed_mask = target_changed_mask[:target] + [True] * len(clip) + target_changed_mask[target:]
                target_junction_mask = target_junction_mask[:target] + [False] * len(clip) + target_junction_mask[
                    target:]
                mark_junction(target_junction_mask, target)
                mark_junction(target_junction_mask, target + len(clip) - 1)
        else:
            logger.warning(f'Unknown OP_TYPE: {rec.info["OP_TYPE"]}')

    if len(changed_mask) == 0:
        logger.warning(f'No changed intervals found for {svid}-{sv_type}')
        return []

    MERGE_TOLERANCE = 5

    queries.extend(
        get_changed_subsequences(new_sequence, changed_mask, junction_mask, MERGE_TOLERANCE, offset, buffer, svid,
                                 chrom, sv_type))
    if not merged_target_sequence:
        queries.extend(
            get_changed_subsequences(new_target_sequence, target_changed_mask, target_junction_mask, MERGE_TOLERANCE,
                                     target_offset, buffer,
                                     svid, chrom, sv_type))

    return queries


def get_changed_subsequences(new_sequence, changed_mask, junction_mask, tolerance, offset, buffer, svid, chrom,
                             sv_type) -> List[QueryReconSubsequence]:
    changed_intervals = []
    queries = []
    # Length of the FULL resulting allele these subsequences are carved from -- a
    # read must be at least this long to contain the whole variant (read mode).
    # Sum the element lengths rather than len(new_sequence): deletions leave empty
    # '' placeholders (and inserted clips can be multi-char), so element count
    # over-estimates the true bp length for deletion-containing alleles.
    result_len = sum(len(s) for s in new_sequence)

    current_index = 0
    for value, group in groupby(changed_mask):
        group_len = len(list(group))
        if value:
            interval_start = current_index
            prev_interval = changed_intervals[-1] if changed_intervals else None

            if prev_interval and current_index - prev_interval[1] <= tolerance:
                prev_interval = changed_intervals.pop()
                interval_start = prev_interval[0]

            changed_intervals.append((interval_start, current_index + group_len))
        current_index += group_len

    for start, stop in changed_intervals:
        adjusted_start = max(0, start - buffer)
        adjusted_stop = min(stop + buffer + 1, len(new_sequence))

        seq_slice = new_sequence[adjusted_start:adjusted_stop]
        junc_slice = junction_mask[adjusted_start:adjusted_stop]

        sequence = ''
        junctions = []
        current_str_idx = 0

        for seq_char, is_junc in zip(seq_slice, junc_slice):
            if is_junc:
                junctions.append(current_str_idx)
            sequence += seq_char
            current_str_idx += len(seq_char)

        junctions = sorted(list(set(junctions)))

        query = QueryReconSubsequence(
            chrom=chrom,
            svtype=sv_type,
            svid=svid,
            sequence=sequence,
            location=adjusted_start + offset,
            length=len(sequence),
            junctions=junctions,
            result_len=result_len,
        )
        queries.append(query)

    return queries
