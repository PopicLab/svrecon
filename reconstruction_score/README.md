# Reconstruction Scorer

**Reconstruction evaluation procedure**

Given a set of called (and stitched) SVs on a genome, we take the fully assembled genomes and measure 
how frequently the called SVs generate actual subsequences in the sample genome.

**Subsequence reconstruction**

Complex SVs may include multiple operations, so we extract the subsequence surrounding the 
positions and targets of the included operations, and implement the resulting subsequence after transformation. We extract the sequence plus and minus a buffer of bps around the endpoints. This results in a subsequence that, if the SV is correctly called, should appear in the indicated chromosome of the sample.

- To search the sample for the subsequence, we load the sample into a minimap2 aligner.
    - Current implementation uses map-pb long read preset for mappy, but this may not be the best setting. The disconnect is that we are seeking the best full match to the query subsequence created by the SV, but aligner objectives can highly score partial matches. We don’t want partial matches, since they will be good partial matches with or without the SV transformations.
    - After finding alignments, we compute the full edit distance between the query and the match, and if the ratio of edit distance to the query length is below a threshold, it’s considered a match, and the SV is considered **correct**.
- We run each SV alone, which means the effects of other SVs in the genome (called or not) are ignored. Thus, the raw location of the SV may not be accurate, so we are searching the full chromosome for occurrences of the resulting sequence.
- This evaluation is usually correct, but it can be wrong in a few circumstances
    - The resulting sequence is in the sample by coincidence, not because of an SV. This likely happens in highly repetitive regions.
    - The SV being measured is near other SV. This happens often with nearby deletions. The buffer regions before and after the SV in question would be wrong if they do not also consider the changes caused by nearby SVs. This only occurs with very close SVs, but they do happen.
- For dispersions, we search for the resulting subsequence at the source and target of the dispersion (if there is also a change at the source).

## Instructions for running:

Run with the following parameters.
- `--reference`: Path to .fa file containing the reference genome
- `--sample`: Path to .fa file containing the sample genome (on which SVs have been called)
- `--calls`: Path to .vcf file containing called SVs in InsilicoSV format
- `--buffer`:  Subsequence context buffer size, number of bps before and after reconstructed SV to compare
- `--gap_file`: (Optional) Tab-delimited file containing regions to omit (e.g., centromere and telomere)

Optionally include the following parameters to generate IGV session xmls to visualize the predictions
- `--config`: (Optional) Groovi call config used to infer other params (will attempt to infer `bam`, `classified`, `reference`, `sample`, and `calls`, so those can be omitted if they are in the config file)
- `--bam`: (Optional) BAM file for generating IGV config
- `--classified`: (Optional) VCF file of groovi-style classified breakpoints for IGV config

Path resolution when inferring from a config file is somewhat sensitive to how we organize files internally, so it may not work in other environments. It's just a shortcut for manually inputting each parameter, so it should be possible to work around.

## Example calls 
- `python score_alignments.py --reference ./data/genome.chr21.fa --sample ./sim_data/sim.combined.fa --calls ./sim_data/sim.vcf --buffer 500`
- `python score_alignments.py --reference /data/refs/refdata-hg19-2.1.0/fasta/genome.fa --sample /data/refs/HG002/hg002v1.1.fasta --calls /data/bert/groovi/vcf_export_debug/results/groovi.vcf --gap_file /data/bert/datasets/hg19.gap.txt`


