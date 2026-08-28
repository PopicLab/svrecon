# svrecon

**svrecon** — an SV callset validation framework based on reconstruction scoring. svrecon validates a callset (`--calls`) by reconstructing each alt sequence, then scoring it against provided long reads (`--bam`), a sample assembly (`--sample`), or both.

## Installation

svrecon is a Python package (Python ≥ 3.11). Install from a clone:

```bash
# pip install
pip install -e .                                    
# or straight from the source repo:
pip install git+https://github.com/PopicLab/svrecon.git
```

This puts an `svrecon` command on your `PATH` (equivalently `python -m svrecon`).

## Tests

A small suite over 20 SVs that insilicoSV simulates in a 200 kb synthetic genome — two instances of 10 classes each (20 total), each checked in both directions: the rebuilt allele must validate, the
unrearranged reference must not.

```bash
pytest
```

The fixture in `tests/data/` is committed. To regenerate it (needs `insilicosv` on `PATH`):

```bash
tests/generators/generate_data.sh
```

## Running

Three things are always required: `--reference`, `--calls`, and at least one validation source.
**The validation mode is inferred from which sources you give**.

```bash
# reads only -- validate each allele against the long reads in the BAM
svrecon --reference ref.fa --calls calls.vcf --bam sample.bam

# assembly only -- validate against the sample assembly
svrecon --reference ref.fa --calls calls.vcf --sample sample.fa

# both -- reads first (Tier 1), assembly as fallback (Tier 2) for alleles no read spans
svrecon --reference ref.fa --calls calls.vcf --bam sample.bam --sample sample.fa
```

Give neither `--bam` nor `--sample` and every call comes back `inconclusive` — nothing was
validated.

For a repeatable run, put the same keys in a YAML config and pass `--config`. **Logs and reports
are written to the config file's directory**, so each experiment is self-contained; without
`--config` they go to the working directory.

```yaml
# experiments/hg002/config.yaml
reference: /data/refs/hg38.fa
calls:     /data/groovi/groovi.vcf
bam:       /data/HG002/hifi.bam
sample:    /data/HG002/hg002v1.1.fasta
report:    json
```

```bash
svrecon --config experiments/hg002/config.yaml   # -> experiments/hg002/svrecon.log, svrecon.report.jsonl
```

Config keys mirror the flag names with underscores (`read_error_threshold`, not
`--read-error-threshold`). Precedence is **CLI flag > config value > default**; setting the same
key in both places is an error.

### All parameters

Inputs:

| Flag | Default | Description |
| --- | --- | --- |
| `--reference` | — | Reference genome `.fa`. Required. |
| `--calls` | — | VCF of called SVs, in InsilicoSV format. Required. |
| `--bam` | — | Long reads aligned to `--reference`. Enables read validation. Indexed on first use if no `.bai`/`.csi` exists. |
| `--sample` | — | Sample genome `.fa`. Enables assembly validation. |
| `--gap-file` | — | Tab-delimited regions to omit (e.g. centromere, telomere). |
| `--config` | — | YAML config; also sets the output directory. |

Validation thresholds:

| Flag | Default | Description |
| --- | --- | --- |
| `--buffer` | `500` | Reference context, in bp, kept on each side of the rebuilt allele. `auto` sizes it per SV to `max(50, 10% of that SV's longest segment)`. |
| `--match-error-threshold` | `0.1` | Max error rate for an assembly / edlib / reference alignment to validate a query. |
| `--read-error-threshold` | `0.1` | Max error rate for a read to validate a query. Keep ≤ `--match-error-threshold` for consistent hit/miss calls. |
| `--location-tolerance` | unbounded | Max bp between a mappy hit's reference start and the expected locus. Bounding it is only safe for well-scaffolded assemblies — per-contig assemblies report contig-local coordinates, so a bound rejects valid alignments. In YAML write infinity as `.inf` (bare `inf` parses as a string). |
| `--max-reads-per-site` | `1000` | Cap on candidate reads gathered per locus; bounds work on deep pileups. |
| `--assembly-forward-match-only` | off | Drop reverse-strand assembly alignments before the CIGAR checks. |
| `--check-reference` | off | Also align each passing allele to the reference; if it validates there too, downgrade the call to `inconclusive`. See below. |

Output:

| Flag | Default | Description |
| --- | --- | --- |
| `--report` | `none` | `json` writes `svrecon.report.jsonl` beside the log. See below. |
| `--verbose` | off | Append the full JSON summary to each per-SV log line. |
| `--n-threads` | half the CPUs | Worker threads for scoring SVs. |
| `--chrom-cache` | temp dir | Where to keep per-chromosome `.mmi` indices between runs. |
| `--plot-first-n` | `0` | Write k-mer dot plots for the first N calls of each SV type, into `<output_dir>/sv_recon_img/`. |
| `--plot-substitute-bases` | off | Replace non-ACGT bases with random ACGT so a plot can be drawn (otherwise that plot is skipped). Plotting only — scoring never sees substituted bases. |
| `--plot-aspect` | `auto` | `auto` keeps the plot square; `equal` is true-to-scale but squeezes asymmetric SVs into a sliver. |
| `--plot-title-svid`, `--plot-title-location`, `--plot-axis-length` | off | Extra labels on the plots. |
| `--classified`, `--igv-prefix` | — | groovi-style classified-breakpoint VCF, and a path prefix, for the generated IGV session XMLs. |

Config-file only (no CLI flag):

| Key | Default | Description |
| --- | --- | --- |
| `max_indel_size` | unset | If set, fail any alignment carrying an indel at least this long. |
| `junction_radius` | unset | If set, additionally check the error rate within this radius of every breakpoint. |
| `min_mappable_fraction` | `0.95` | A query with less than this fraction of mappable (ACGT) bases is inconclusive rather than pass/fail. Set to `null` to disable. Only bites above `1 - error threshold` — see below. |
| `auto_buffer_min`, `auto_buffer_fraction` | `50`, `0.1` | Floor and fraction used by `buffer: auto`. |
| `min_edlib_query` | `5000` | Queries at least this long skip the edlib fallback (cost bound). |
| `edlib_fallback_max_tolerance` | `10000000` | Search-window radius cap for the edlib fallback. A cost bound on an otherwise unbounded O(n·m) search, kept separate from `location_tolerance`. |

## Procedure

For each SV, svrecon reconstructs the sequence the call implies and asks whether that sequence
exists in the sample.

**1. Reconstruct.** A complex SV may carry several operations, so each SV yields one or more
**queries**: the reference spanning the operations' positions and targets, transformed as the call
describes, plus `--buffer` bp of reference context on each side. If the call is correct, that
query should appear in the sample. Dispersions produce a query at both the source and the target.
Each SV is reconstructed in isolation, so nearby SVs — called or not — are ignored.

The reconstruction process parses a VCF in insilicoSV format to produce an initial set of disjoint
reference segments, as well as a sequence of deletion, insertion, and inversion operations, where a
complex record can produce more than one. To apply them unambiguously, these operations operations are ordered:

1. Rightmost-first, so an applied edit never shifts a still-pending operation's coordinates.
2. At a tied position, invert/delete before insert -- they need an untouched segment boundary,
   which a same-position insert would shift.
3. Among tied inserts, by descending INSORD (insilicoSV's insertion-order field) -- each insert at a
   shared position lands to the left of ones already placed there, thus processing operations in reverse results in placements of increasing INSORD order, left to right.

**2. Align.** Each query is aligned against the configured sources in tier order, stopping at the
first pass:

- **reads** (`--bam`) — fetch reads overlapping the locus, keep those at least as long as the
  query, and align the query into each with edlib, trying both orientations. A read can span the
  reference source and still contain only part of the rebuilt allele, so covering the source
  region is *not* sufficient; the read molecule must be at least the full allele length.
- **assembly** (`--sample`) — minimap2/mappy (`map-hifi`, `--eqx`) against a per-chromosome index,
  then filter hits by `--location-tolerance`.
- **edlib fallback** (`--sample`) — for queries shorter than `min_edlib_query`, an
  expanding-window edlib search of the sample bytes, catching short alleles mappy missed.

**3. Check the CIGAR.** An alignment alone is not enough — every configured check must pass on it.
This catches alignments whose overall error is diluted below threshold by the flanking buffer but
that are wrong precisely where the SV rearranges the sequence.

- *sequence similarity* (always) — total error over the whole query ≤ the threshold.
- *no large indels* (if `max_indel_size` is set) — no single indel that long.
- *junctions* (if `junction_radius` is set) — local error within that radius of every breakpoint.
- *mappable* (on by default, `min_mappable_fraction: 0.95`) — proportion of query sequence that must be mappable bases (ACGT) for the query to not be marked inconclusive


A check returns `pass`, `fail`, or `inconclusive`, and the alignment's verdict is the roll-up:
`inconclusive` if any check could not judge the query, else `fail` if any failed, else `pass`.

**4. Classify.** Each query gets a status and a reason:

| status | reason | meaning |
| --- | --- | --- |
| `pass` | `pass` | Aligned, and every check passed. |
| `fail` | `cigar_failed` | **Aligned**, then a check rejected it — the sample contradicts the call. |
| `fail` | `other` | Nothing aligned within the error budget. |
| `inconclusive` | `other` | Untestable: no read long enough to span the allele, a check could not judge the query (too little of it mappable), the allele also matches the reference, or no validation source was configured. |

`cigar_failed` is exactly the marker that the query **aligned**; the report exposes it directly as
`"aligned": true`. When several scorers run, the decisive result is the most settling one:
`pass` > aligned `fail` > `inconclusive` > unaligned `fail`, ties broken by lower error.

**5. Roll up to the SV.** `hit` if every query passed, `miss` if any query failed, else
`inconclusive`. Precision counts only conclusive calls — `hits / (hits + misses)` — so
inconclusive calls are excluded from the denominator and reported in their own column.

Hits also carry an **evaluation tier**: `read` (Tier 1 — a single real read contained the allele;
the strongest evidence, independent of assembly quality) or `assembly` (Tier 2 — confirmed only
against the reconstruction, via mappy or edlib). Even a T2T-grade assembly is itself a
reconstruction, which is why a read-confirmed call ranks higher. The score table's `read_hits` and
`assembly_hits` columns are the per-tier counts.

Two failure modes this cannot rule out on its own: an allele present in the sample by coincidence
rather than by an SV (likely in repetitive regions — see `--check-reference`), and a call adjacent
to another SV, which corrupts its buffer.

## Reference-ambiguity check (`--check-reference`)

Validation confirms the allele exists in the sample — but in repetitive or segmental-duplication
regions the allele can exist in the **reference** too, in which case finding it in the sample says
nothing about whether the SV occurred. With `--check-reference`, each *passing* allele is also
aligned to the reference; if it validates there as well, the call is downgraded from a hit to
`inconclusive`, flagged `"reference_ambiguous": true` in the report.

The check uses mappy, not edlib, because a large or compound allele needs chained alignment —
plain edlib cannot align a multi-junction allele and would spuriously report "absent". It
therefore builds a second per-chromosome aligner set over the reference (extra build time and
memory), so use it deliberately (auditing suspicious large calls) rather than on every run. It
caught, for example, a 29 kb chr16 dupINVdup whose allele maps to the reference at 0.019 error —
below its 0.034 assembly match — i.e. a coincidental, non-SV-specific hit.

## Per-SV report (`--report json`)

The log's per-SV lines give one outcome per call, plus either the tier (on a hit) or a
`status:reason` tag per query (on a miss or inconclusive). For deeper inspection — *which checks
ran, which passed, and what each measured* — pass `--report json`. It writes one JSON record per
SV to `svrecon.report.jsonl` beside the log. This is additive: the log and score table are
identical whether or not it is enabled.

Each record holds only what the *evaluation* concluded; join back to the call VCF on `svid` for
coordinates, types, and operations.

```jsonc
{"svid": "sv10319",
 "svtype": "dupINVdup",          // redundant with the VCF; included for quick scanning
 "outcome": "hit",               // hit | miss | inconclusive
 "tier": "assembly",             // read | assembly | miss
 "query_validations": [          // one entry per reconstructed query
  {"chrom": "chr3",              // chromosome the query maps to (pairs with ref_start)
   "ref_start": 187135426,       // reference coord of the query window's left anchor
   "status": "pass",             // decisive verdict: pass | fail | inconclusive
   "reason": "pass",             // pass | cigar_failed | other
   "aligned": true,              // whether an alignment was found at all
   "reference_ambiguous": false, // true if the allele also validated against the reference
   "source": "assembly",         // reads | assembly | edlib | null (which scorer decided)
   "strand": 1,                  // 1 | -1 | null (null for reads -- randomly oriented)
   "validations": [              // every scorer that ran, in order; the decisive one is ranked, not last
    {"source": "assembly", "passed": true, "status": "pass", "reason": "pass", "aligned": true,
     "cigar": "500=3X2I495=",      // the decisive alignment's CIGAR, for debugging
     "lowest_pass_error": 0.0696,  // error of the adopted alignment (1.0 if none passed)
     "lowest_error": 0.0696,       // best error seen, pass or fail
     "strand": 1,
     "checks": [                   // one per configured CIGAR check, in configured order
      {"status": "pass", "detail": "overall error 0.0696 <= 0.1"},
      {"status": "pass", "detail": "junction errors {0:0.0, 812:0.01}, worst 812:0.01 <= 0.1"}]}],
   "ambiguity_validations": []}]}  // same shape, against the reference (--check-reference only)
```

`ref_start` is the window anchor (≈ breakpoint − buffer), not an exact breakpoint, and a
multi-operation complex SV produces one `query_validations` entry per reconstructed query, not per
VCF record.

## Benchmark

[`workflows/svrecon_benchmark.ipynb`](workflows/svrecon_benchmark.ipynb) measures what the CIGAR
checks buy you. It simulates 120 SVs (8 in-place types x 3 size classes x 5) on hg38 chr21 with
insilicoSV, then scores two callsets against one assembly:

| arm | callset | assembly | ideal |
| --- | --- | --- | --- |
| positive | seed 0 | seed 0 | 120/120 hits |
| negative | seed 1 | seed 0 | 0/120 hits |

Each arm is scored under three check configurations — `similarity-only`,
`similarity-no-large-indels` (`max_indel_size: 50`), and `similarity-junction`
(`junction_radius: 100`) — so a hit in the negative arm is a false positive attributable to that
configuration. All three reach 120/120 on the positive arm; on the negative arm they leave 12, 12,
and 0 false positives respectively. Run it from the repo root with `insilicosv` and `svrecon` installed; the configs
live in `workflows/`, and the last cell deletes everything generated.

## Notes

- Outputs go to the `--config` file's directory; without a `--config`, to the working directory.
- Per-chromosome `.mmi` indices are cached (in temp space unless `--chrom-cache` is set) and
  rebuilt automatically when the source FASTA changes.
- Dot plots use wotplot, which accepts only `A`/`C`/`G`/`T`; a query containing `N` is skipped
  unless `--plot-substitute-bases` is given.
