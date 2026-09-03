#!/usr/bin/env bash

# Script to reproducibly generate BAM and FASTA files for test.
set -euo pipefail
cd "$(dirname "$0")/../data"

# Clear stale outputs. The .fai/.bai matter: pysam trusts a stale index and silently
# reads the wrong length.
rm -f assembly.* reference.fa* reads.bam*

# Run InsilicoSV
PYTHONPATH=.. python ../generators/generate_small_genome.py
insilicosv -c insilicoSV.yaml

# Clean up InsilicoSV outputs, rename sim* to assembly*
for f in sim.*; do mv "$f" "assembly.${f#sim.}"; done
mv assembly.hapA.fa assembly.fa 
rm -f assembly.hapB.fa assembly.divergence.fa assembly.novel_insertions.fa assembly.stats.txt

PYTHONPATH=.. python ../generators/generate_small_bam.py
