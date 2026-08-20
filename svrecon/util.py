"""Small shared helpers: sequence ops, FASTA/VCF loading, IGV session export."""
import os
from collections import defaultdict
from typing import Dict, List, Optional, Tuple, Union

import pysam
from intervaltree import IntervalTree
from pysam import VariantRecord

_COMPLEMENT_TRANS = str.maketrans('ACGTNacgtn', 'TGCANtgcan')

def clamp(x: float, lo: float, hi: float) -> float:
    """Clamps x to the closed interval [lo, hi]."""
    return max(lo, min(hi, x))


def merge_intervals(intervals: List[List[int]]) -> List[List[int]]:
    """Merges overlapping [start, end] intervals into their union."""
    merged = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return merged

def reverse_complement(seq: Union[str, List[str]]) -> str:
    """Reverse complement of a DNA sequence (unknown bases -> N). Accepts a str or
    a list of single-character strings; always returns a str."""
    if not isinstance(seq, str):
        seq = ''.join(seq)
    return seq.translate(_COMPLEMENT_TRANS)[::-1]


def load_fasta_to_bytes(filename: str, chroms) -> Dict[str, bytearray]:
    """Loads target sequences into a dictionary of bytearrays."""
    data = {}
    with pysam.FastaFile(filename) as f:
        fasta_refs = f.references
        for chrom in chroms:
            chrom_lower = chrom.lower()
            for ref in fasta_refs:
                ref_lower = ref.lower()
                if ref_lower == chrom_lower or ref_lower.startswith(f"{chrom_lower}_"):
                    sequence_string = f.fetch(reference=ref)
                    data[ref] = bytearray(sequence_string, 'ascii')
    return data


def get_start_stop(rec: VariantRecord) -> Tuple[int, int]:
    """Obtains VCF start/stop coordinates. If SVLEN is available, calculate stop (due to pysam bug in shifting stop)"""
    start = rec.start
    stop = start + rec.info['SVLEN'] if rec.info.get('SVLEN') else rec.stop
    return start, stop


def load_exclude_list(gap_file: str) -> Dict[str, IntervalTree]:
    exclude_list = defaultdict(IntervalTree)
    with open(gap_file, 'r') as f:
        for line in f:
            row = line.strip().split()
            chrom = row[1]
            start, stop = int(row[2]), int(row[3])
            region_type = row[7]
            exclude_list[chrom][start:stop] = region_type
    return exclude_list


def group_variants_by_id(vcf_path: str, gap_file: Optional[str] = None) -> Dict[str, List[VariantRecord]]:
    """Group a callset VCF's records by SVID, optionally dropping any SV with a record
    (its own span, or its TARGET) overlapping an excluded region (e.g. centromere/telomere)."""
    grouped_variants: Dict[str, List[VariantRecord]] = defaultdict(list)
    for rec in pysam.VariantFile(vcf_path).fetch():
        svid = rec.info.get('SVID')
        if svid:
            grouped_variants[svid].append(rec)

    if not gap_file:
        return grouped_variants

    exclude_list = load_exclude_list(gap_file)
    filtered_variants: Dict[str, List[VariantRecord]] = {}
    for svid, records in grouped_variants.items():
        allowed = True
        for rec in records:
            target_chrom = rec.info['TARGET_CHROM'] if 'TARGET_CHROM' in rec.info else None
            if exclude_list[rec.chrom].overlap(rec.start, rec.stop) or \
                    target_chrom and exclude_list[target_chrom].overlap(rec.info['TARGET'],
                                                                        rec.info['TARGET'] + 1):
                allowed = False
                break
        if allowed:
            filtered_variants[svid] = records

    return filtered_variants


def export_igv_session(calls, bam, classified, timestamp, igv_prefix):
    import xml.etree.ElementTree as ET
    output_filename = f'./igv_sessions/session_{timestamp}.xml'

    session = ET.Element('Session', genome='hg19', version='8')
    resources = ET.SubElement(session, 'Resources')

    if calls:
        ET.SubElement(resources, 'Resource', path=igv_prefix + calls, type='vcf')
    if classified:
        ET.SubElement(resources, 'Resource', path=igv_prefix + classified, type='vcf')
    if bam:
        bam_path = igv_prefix + bam
        ET.SubElement(resources, 'Resource', path=bam_path, type='bam')
        bam_panel = ET.SubElement(session, 'Panel', name='Alignments', height='600')
        align_track = ET.SubElement(bam_panel, 'Track', clazz='org.broad.igv.sam.AlignmentTrack',
                                    displayMode='EXPANDED', id=bam_path, name=os.path.basename(bam_path),
                                    visible='true')
        ET.SubElement(align_track, 'RenderOptions', colorOption='READ_STRAND', duplicatesOption='FILTER',
                      groupByOption='LINKED', hideSmallIndels='true', linkByTag='READNAME',
                      linkedReads='true', smallIndelThreshold='2')

    os.makedirs(os.path.dirname(output_filename), exist_ok=True)
    tree = ET.ElementTree(session)
    tree.write(output_filename, encoding='utf-8', xml_declaration=True)
    print(f'File saved to {output_filename}')
