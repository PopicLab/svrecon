"""Generate the synthetic reference genome that insilicoSV simulates against."""
import random
from pathlib import Path

import pysam

CHROM = 'chrT'
CHROM_LEN = 200_000  # room for 20 SVs at min_intersv_dist 5000
SEED = 0
FASTA_LINE_WIDTH = 60
DATA_DIR = Path(__file__).parents[1] / 'data'
REFERENCE_FASTA = DATA_DIR / 'reference.fa'


def random_sequence(length: int, seed: int) -> str:
    rng = random.Random(seed)
    return ''.join(rng.choice('ACGT') for _ in range(length))


def write_genome(path=REFERENCE_FASTA) -> str:
    """Write the reference FASTA and its .fai. Returns the path."""
    sequence = random_sequence(CHROM_LEN, SEED)
    with open(path, 'w') as fasta:
        fasta.write(f'>{CHROM}\n')
        for offset in range(0, len(sequence), FASTA_LINE_WIDTH):
            fasta.write(sequence[offset:offset + FASTA_LINE_WIDTH] + '\n')
    pysam.faidx(str(path))
    return str(path)


if __name__ == '__main__':
    print(f'wrote {write_genome()} ({CHROM}, {CHROM_LEN} bp)')
