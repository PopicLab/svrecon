"""Small shared helpers: sequence ops, FASTA loading, IGV session export."""
import os
from typing import Dict, List, Tuple, Union

import pysam
from pysam import VariantRecord


_RC_TRANS = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


def reverse_complement(seq: Union[str, List[str]]) -> str:
    """Reverse complement of a DNA sequence (unknown bases -> N). Accepts a str or
    a list of single-character strings; always returns a str."""
    if not isinstance(seq, str):
        seq = ''.join(seq)
    return seq.translate(_RC_TRANS)[::-1]


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
    """Converts VCF start/stop coordinates to Python style"""
    start = rec.start - 1
    stop = rec.stop
    return start, stop


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
