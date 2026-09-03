"""Generate one supporting read per simulated SV, cut from the assembly insilicoSV produced."""
from dataclasses import dataclass
from typing import List

import pysam

from generate_small_genome import CHROM, CHROM_LEN, DATA_DIR, REFERENCE_FASTA
from svrecon.utils import get_start_stop, load_grouped_variants_from_vcf

ASSEMBLY_VCF = DATA_DIR / 'assembly.vcf'
ASSEMBLY_FASTA = DATA_DIR / 'assembly.fa'
ASSEMBLY_CHROM = f'{CHROM}_hapA'
READS_BAM = DATA_DIR / 'reads.bam'

# Read info
READ_BUFFER = 1_000
MAPPING_QUALITY = 60
HEADER = pysam.AlignmentHeader.from_dict({'HD': {'VN': '1.6', 'SO': 'coordinate'},
                                          'SQ': [{'SN': CHROM, 'LN': CHROM_LEN}]})


def load_sequence(fasta_path, chrom: str) -> str:
    with pysam.FastaFile(str(fasta_path)) as fasta:
        return fasta.fetch(chrom).upper()


# --- read the simulated SVs off the VCF ---

@dataclass
class SimulatedSV:
    """One SV from assembly.vcf, and the assembly sequence replacing its reference span."""
    svid: str
    svtype: str
    ref_start: int
    ref_end: int
    alt_sequence: str


def locate_in_assembly(sequence: str, assembly: str, search_from: int = 0) -> int:
    """Index of ``sequence`` in the assembly."""
    position = assembly.find(sequence, search_from)
    assert position >= 0, 'sequence is missing from the assembly'
    assert assembly.find(sequence, position + 1) < 0, 'sequence is not unique in the assembly'
    return position


def load_simulated_svs(vcf_path, reference: str, assembly: str) -> List[SimulatedSV]:
    """Locate reconstructed SVs by their surrounding (flanking) regions"""
    simulated_svs = []
    for svid, records in load_grouped_variants_from_vcf(str(vcf_path)).items():
        spans = [get_start_stop(rec) for rec in records]
        ref_start = min(start for start, _ in spans)
        ref_end = max(stop for _, stop in spans)

        left_flank_at = locate_in_assembly(reference[ref_start - READ_BUFFER:ref_start], assembly)
        right_flank_at = locate_in_assembly(reference[ref_end:ref_end + READ_BUFFER], assembly,
                                            left_flank_at)
        simulated_svs.append(SimulatedSV(
            svid=svid, svtype=records[0].info['SVTYPE'], ref_start=ref_start, ref_end=ref_end,
            alt_sequence=assembly[left_flank_at + READ_BUFFER:right_flank_at]))
    return sorted(simulated_svs, key=lambda sv: sv.ref_start)


# --- generate reads ---

def make_read(name: str, ref_start: int, sequence: str) -> pysam.AlignedSegment:
    """Mapped full-length match spanning ``[ref_start, ref_start + len(sequence))``.
    The sequence is stored verbatim, so it need not match the reference the CIGAR claims."""
    read = pysam.AlignedSegment(HEADER)
    read.query_name = name
    read.reference_id = 0
    read.reference_start = ref_start
    read.mapping_quality = MAPPING_QUALITY
    read.query_sequence = sequence
    read.cigartuples = [(pysam.CMATCH, len(sequence))]
    return read


def make_sv_read(sv: SimulatedSV, reference: str) -> pysam.AlignedSegment:
    """The read supporting ``sv``: its assembly sequence, READ_BUFFER bp of reference each side."""
    name = f'{sv.svtype}_read'
    read_start = sv.ref_start - READ_BUFFER
    left_flank = reference[read_start:sv.ref_start]
    right_flank = reference[sv.ref_end:sv.ref_end + READ_BUFFER]
    sequence = left_flank + sv.alt_sequence + right_flank
    return make_read(name, read_start, sequence)


# --- write the bam ---

def write_bam(path, reads: List[pysam.AlignedSegment]) -> str:
    """Coordinate-sort, write, index. Returns the path; a .bai sits beside it."""
    with pysam.AlignmentFile(str(path), 'wb', header=HEADER) as bam:
        for read in sorted(reads, key=lambda read: read.reference_start):
            bam.write(read)
    pysam.index(str(path))
    return str(path)


if __name__ == '__main__':
    reference = load_sequence(REFERENCE_FASTA, CHROM)
    assembly = load_sequence(ASSEMBLY_FASTA, ASSEMBLY_CHROM)
    simulated_svs = load_simulated_svs(ASSEMBLY_VCF, reference, assembly)
    reads = [make_sv_read(sv, reference) for sv in simulated_svs]

    for sv, read in zip(simulated_svs, reads):
        print(f'  {sv.svtype:10} {sv.svid:5} source [{sv.ref_start},{sv.ref_end}) '
              f'alt {len(sv.alt_sequence):5} bp  read {read.query_length} bp '
              f'@ {read.reference_start}')
    print(f'wrote {write_bam(READS_BAM, reads)} ({len(reads)} reads on {CHROM})')
