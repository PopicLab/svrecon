# svrecon

**svrecon** — structural-variant reconstruction scoring. Validate called SVs by rebuilding each variant's alt allele and checking whether real sequence supports it, against a sample assembly (`assembly`), long reads (`reads`), or reads-first with assembly fallback (`both`).

## Instructions for running:

Run with the following parameters.
- `--reference`: Path to .fa file containing the reference genome
- `--sample`: Path to .fa file containing the sample genome (on which SVs have been called)
- `--calls`: Path to .vcf file containing called SVs in InsilicoSV format
- `--buffer`:  Subsequence context buffer size, number of bps before and after reconstructed SV to compare
- `--gap_file`: (Optional) Tab-delimited file containing regions to omit (e.g., centromere and telomere)
- `--location_tolerance`: (Optional) Max bp between a mappy hit's reference start and the expected SV locus for the hit to count; default **unbounded**. Bounding it is only safe for well-scaffolded assemblies — for per-contig/unscaffolded assemblies, whose hit coordinates are contig-local rather than genomic, a bounded tolerance rejects valid matches.
- `--chrom_cache`: (Optional) Directory to save/load the per-chromosome `.mmi` indices (defaults to temporary space).
- `--report`: (Optional) `none` (default) or `json`. `json` writes a per-SV evaluation sidecar (see "Per-SV evaluation report" below) alongside the log. Off by default; the log and score table are identical either way.

Read-based evaluation parameters (see "Read-based evaluation" below):
- `--eval_mode`: `assembly` (default), `reads`, or `both`. `reads` validates each call against the long reads in `--bam` instead of the assembly; `both` is **reads-first** — it consults the reads (Tier 1) and only falls back to the assembly (Tier 2) for events no single read can span.
- `--bam`: BAM of long reads aligned to `--reference` (required for `reads`/`both`).
- `--read_error_threshold`: max edlib error for a read to confirm a reconstruction (default 0.1; keep ≤ the assembly threshold).
- `--min_read_support`: min number of spanning reads that must clear the threshold (default 1).
- `--max_reads_per_site`: cap on candidate reads gathered per locus (default 1000).

Optionally include the following parameters to generate IGV session xmls to visualize the predictions
- `--config`: (Optional) Groovi call config used to infer other params (will attempt to infer `bam`, `classified`, `reference`, `sample`, and `calls`, so those can be omitted if they are in the config file). See notes.
- `--bam`: (Optional) BAM file for generating IGV config
- `--classified`: (Optional) VCF file of groovi-style classified breakpoints for IGV config
- `--igv_prefix`: (Optional) Path prefix prepended to file paths in the generated IGV session XMLs (e.g., a local mount point for files that live on a remote server).

Path resolution when inferring from a config file is somewhat sensitive to how we organize files internally, so it may not work in other environments. It's just a shortcut for manually inputting each parameter, so it should be possible to work around.

## Example calls 
- `python score_alignments.py --reference ./data/genome.chr21.fa --sample ./sim_data/sim.combined.fa --calls ./sim_data/sim.vcf --buffer 500`
- `python score_alignments.py --reference /data/refs/refdata-hg19-2.1.0/fasta/genome.fa --sample /data/refs/HG002/hg002v1.1.fasta --calls /data/bert/groovi/vcf_export_debug/results/groovi.vcf --gap_file /data/bert/datasets/hg19.gap.txt`




**Reconstruction evaluation procedure**

Given a set of called (and stitched) SVs on a genome, we take the fully assembled genomes and measure 
how frequently the called SVs generate actual subsequences in the sample genome.

**Subsequence reconstruction**

Complex SVs may include multiple operations, so we extract the subsequence surrounding the 
positions and targets of the included operations, and implement the resulting subsequence after transformation. We extract the sequence plus and minus a buffer of bps around the endpoints. This results in a subsequence that, if the SV is correctly called, should appear in the indicated chromosome of the sample.

- To search the sample for the subsequence, we load the sample into a minimap2 aligner.
    - Current implementation uses the map-hifi preset for mappy, but this may not be the best setting. The disconnect is that we are seeking the best full match to the query subsequence created by the SV, but aligner objectives can highly score partial matches. We don’t want partial matches, since they will be good partial matches with or without the SV transformations.
    - After finding alignments, we compute the full edit distance between the query and the match, and a match requires **both** (a) the ratio of edit distance to the query length below a threshold, **and** (b) the SV's breakpoint junctions to validate — the reconstructed adjacencies must align cleanly within a ±300 bp window (`JUNCTION_VALIDATION_WINDOW`) around each junction. Only then is the SV considered **correct**. The junction check catches alignments whose overall error is diluted below threshold by the flanking buffer but that are wrong precisely at the novel adjacency.
- We run each SV alone, which means the effects of other SVs in the genome (called or not) are ignored. Thus, the raw location of the SV may not be accurate, so we are searching the full chromosome for occurrences of the resulting sequence.
- This evaluation is usually correct, but it can be wrong in a few circumstances
    - The resulting sequence is in the sample by coincidence, not because of an SV. This likely happens in highly repetitive regions.
    - The SV being measured is near other SV. This happens often with nearby deletions. The buffer regions before and after the SV in question would be wrong if they do not also consider the changes caused by nearby SVs. This only occurs with very close SVs, but they do happen.
- For dispersions, we search for the resulting subsequence at the source and target of the dispersion (if there is also a change at the source).

**Read-based evaluation (`--eval_mode reads` / `both`)**

Assembly-based scoring is only as good as the sample assembly. For samples whose
assembly is fragmented, divergent, or coarsely scaffolded, a correct SV can fail
to reconstruct simply because the assembly is a poor target. Read-based mode
sidesteps the assembly: it validates each reconstructed allele directly against
the long reads in the BAM.

How it works:
For each SV, we build the same reconstructed allele (`simulate_subsequences`),
  then gather candidate reads: any read whose alignment (primary or
  supplementary) overlaps the SV position.

Outcomes. Each SV is one of:
- **hit** — a read contains the full reconstructed allele within the error
  threshold (and its junctions validate).
- **miss** — reads that are long enough to contain the whole resulting allele
  exist, but none match it. The reads genuinely contradict the call.
- **inconclusive** — reads overlap the locus, but none is long enough to
  contain the full resulting allele, so a single read cannot confirm or refute
  it. These are excluded from the precision denominator (precision = hits /
  (hits + misses)), and reported separately.

A read can fully span the reference source (reach both breakpoints) while
containing only part of the resulting allele, so source-region coverage is
not a valid test of whether the read spans the called variant. We therefore require a read's
molecule length ≥ the full resulting-allele length before
it can validate. Events whose resulting allele is longer than any read are
inherently inconclusive under single-read validation.

Miss / inconclusive diagnostics are appended to each non-hit log line as a
per-subsequence list, one tag per reconstructed subsequence:
- `no_reads` — no read overlaps the locus (coverage gap).
- `inconclusive` — reads overlap but none span the full resulting allele.
- `over_error_threshold:<err>` — a spanning read aligned, 1×–2× threshold.
- `over_error_threshold:aborted` — spanning reads aligned worse than 2× threshold (edlib aborted).
- `junction_failed:<err>` — aligned within the overall threshold but the breakpoint window failed.

Hit log lines instead carry the **evaluation tier** and the detailed source:
`read` (Tier 1 — a single real read contained the allele; the strongest evidence,
independent of assembly quality) or `assembly` (Tier 2 — confirmed only against the
reconstructed assembly, via `assembly`/`edlib`). The score table's `read_hits` and
`assembly_hits` columns are the per-tier hit counts, and the summary logs a
`Hits by tier` line. In `assembly` mode nothing is inconclusive and every hit is
Tier 2 — behavior is unchanged.

**Two-tier rationale (`both` mode).** Even a T2T-grade assembly is itself a
reconstruction, so a call confirmed by an actual read is stronger evidence than
one confirmed only against an assembly. `both` mode therefore tries reads first and
attributes each hit to the highest tier that confirmed it, falling back to the
assembly only for alleles longer than any single read. In `both` mode, an SV is
`inconclusive` only when neither tier could test it (no read spans the allele and
the assembly produced no candidate alignment); if either tier aligns a candidate
that disagrees, it is a genuine `miss`.

**Per-SV evaluation report (`--report json`)**

The score table and log summarize the run; the log's per-SV lines give one outcome
(plus tier or a compact diagnostic) per call. For deeper inspection — *which checks
ran, which passed, and the local errors at each breakpoint* — pass `--report json`.
It writes a JSONL sidecar next to the log (`logs/<logbasename>.eval.jsonl`), one
record per SV. This is additive: the log and score table are byte-for-byte identical
whether or not it is enabled.

Each record holds only what the *evaluation* concluded; join back to the call VCF on
`svid` for coordinates, types, and operations (which are not duplicated here):

```json
{"svid": "sv10319",
 "svtype": "dupINVdup",         // redundant with the VCF; included for quick scanning
 "outcome": "hit",              // hit | miss | inconclusive
 "tier": "assembly",            // read | assembly | null (null unless hit)
 "segments": [                  // one entry per reconstructed subsequence
   {"chrom": "chr3",            // chromosome the subsequence maps to (pairs with ref_start)
    "ref_start": 187135426,     // reference coord of the subsequence window's left anchor
    "status": "pass",           // pass | fail | inconclusive
    "reason": "pass",           // pass | junction_failed | over_error_threshold |
                                //   no_reads | inconclusive | no_aligner | no_alignment_in_window
    "source": "assembly",       // reads | assembly | edlib | null (tier that decided it)
    "error": 0.0696,            // number, 4 sig figs: winning error if pass, best failing error if
                                //   fail, null if untestable. Tiny rates keep precision and serialize
                                //   in exponent form (e.g. 8.6e-06).
    "junctions": [              // per-junction results in the order checked; truncated at the
      {"error": 0.0, "passed": true}, ...]}]}  // first failure. Empty if no junction in scope.
```

Note the `ref_start` is the window anchor (≈ breakpoint − buffer), not an exact
breakpoint, and a multi-operation complex SV produces one `segments` entry per
reconstructed subsequence, not per VCF record.

**Notes**

- Extracting parameters from groovi config files is very dependent on expected folder structures and likely to break.
  - It expects remote servers to be mounted locally
  - It expects a fasta file to be located in a sibling folder of the .bam file, which typically only occurs for synthetic training data.
  - It expects the call file to be `groovi.vcf` in the config file's results directory
  - Each of these parameters should be overrided if these conditions aren't met. 
- The tool writes to a `logs` directory in the working directory, creating it if absent.
- If the tool doesn't find a `.mmi` index file attached to the `.fa` sample assembly, it creates one as a cache and stores it in temporary space