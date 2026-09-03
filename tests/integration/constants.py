"""Every constant the integration suite and its generators share: paths, names, sizes."""
from pathlib import Path

DATA_DIR = Path(__file__).parent / 'data'

# The synthetic reference genome
CHROM = 'chrT'
CHROM_LEN = 200_000  # room for 20 SVs at min_intersv_dist 5000 # NOTE: increase if test suite expands
REFERENCE_SEED = 0 # seed for random genome simulation
FASTA_LINE_WIDTH = 60
REFERENCE_FASTA = DATA_DIR / 'reference.fa'

# synthetic genome and operations output
ASSEMBLY_FASTA = DATA_DIR / 'assembly.fa'
ASSEMBLY_VCF = DATA_DIR / 'assembly.vcf'
ASSEMBLY_CHROM = f'{CHROM}_hapA'

# read data
READS_BAM = DATA_DIR / 'reads.bam'
READ_BUFFER_SIZE = 1_000
MAPPING_QUALITY = 60

# reconstruction queries
QUERY_BUFFER_SIZE = 500
PLACEHOLDER_GRAMMAR = 'dummy'
