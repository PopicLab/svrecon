"""Generate the synthetic reference genome that insilicoSV simulates against."""
import random

import pysam

from constants import CHROM, CHROM_LEN, FASTA_LINE_WIDTH, REFERENCE_FASTA, REFERENCE_SEED


def generate_random_sequence(length: int, seed: int) -> str:
    rng = random.Random(seed)
    return ''.join(rng.choice('ACGT') for _ in range(length))


def write_genome(path=REFERENCE_FASTA) -> str:
    """Write the reference FASTA and its .fai. Returns the path."""
    sequence = generate_random_sequence(CHROM_LEN, REFERENCE_SEED)
    with open(path, 'w') as fasta:
        fasta.write(f'>{CHROM}\n')
        for offset in range(0, len(sequence), FASTA_LINE_WIDTH):
            fasta.write(sequence[offset:offset + FASTA_LINE_WIDTH] + '\n')
    pysam.faidx(str(path))
    return str(path)


if __name__ == '__main__':
    print(f'wrote {write_genome()} ({CHROM}, {CHROM_LEN} bp)')
