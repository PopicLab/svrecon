# svrecon

**svrecon** — a complex SV callset validation framework based on reconstruction scoring. svrecon validates a callset (`--calls`) by reconstructing each alt sequence from a composable sequence of operations — capable of expressing any complex variant — then scoring it against provided long reads (`--bam`), a sample assembly (`--sample`), or both.

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

```bash
pytest
```

Reconstruction verified against each SV type's grammar, and the CIGAR checks over synthetic
alignments.

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
validated. `--calls` must be a VCF in [insilicoSV](https://github.com/PopicLab/insilicoSV) format.

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
svrecon --config experiments/hg002/config.yaml   # -> experiments/hg002/{svrecon.log,svrecon.summary.json,svrecon.report.jsonl}
```

Config keys mirror the flag names with underscores (`read_error_threshold`, not
`--read-error-threshold`). Precedence is **CLI flag > config value > default**; setting the same
key in both places is an error.

### All parameters

Inputs:

| Flag | Default | Description |
| --- | --- | --- |
| `--reference` | — | Reference genome `.fa`. Required. |
| `--calls` | — | VCF of called SVs, in [insilicoSV](https://github.com/PopicLab/insilicoSV) format. Required. |
| `--bam` | — | Long reads aligned to `--reference`. Enables read validation. Indexed on first use if no `.bai`/`.csi` exists. |
| `--sample` | — | Sample genome `.fa`. Enables assembly validation. |
| `--gap-file` | — | Tab-delimited regions to omit (e.g. centromere, telomere). |
| `--config` | — | YAML config; also sets the output directory. |

Validation thresholds:

| Flag | Default | Description |
| --- | --- | --- |
| `--buffer` | `100` | Reference context, in bp, kept on each side of the rebuilt allele. `auto` sizes it per SV to `max(50, 10% of that SV's longest segment)`. |
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
| `--verbose` | on | Append the full JSON summary to each per-SV log line. |
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
| `max_contiguous_error` | `50` | Minimum contiguous length of a failing insertion, deletion, or soft clip. Set to `null` to disable. |
| `max_window_size`, `max_window_error` | `100`, `25` | If both are set, fail any alignment where some window of `max_window_size` alignment columns holds at least `max_window_error` error columns (mismatches, insertions, deletions, clips). Must be set together; set both to `null` to disable. |
| `min_mappable_fraction` | `0.95` | A query with less than this fraction of mappable (ACGT) bases is inconclusive rather than pass/fail. Set to `null` to disable. Only bites above `1 - error threshold` — see below. |
| `max_unmappable_size` | `50` | A query with a contiguous non-ACGT run at least this long is inconclusive rather than pass/fail. Set to `null` to disable. |
| `auto_buffer_min`, `auto_buffer_fraction` | `50`, `0.1` | Floor and fraction used by `buffer: auto`. |
| `min_edlib_query` | `5000` | Queries at least this long skip the edlib fallback (cost bound). |
| `edlib_fallback_max_tolerance` | `10000000` | Search-window radius cap for the edlib fallback. A cost bound on an otherwise unbounded O(n·m) search, kept separate from `location_tolerance`. |

## Input callset (VCF)

Records are grouped into SVs by their `SVID`; a record without one is a single-record SV of its
own, keyed `simple_{n}` by its position in the callset.

| field | required on | meaning |
| --- | --- | --- |
| `SVID` | every record of a multi-record SV | The key its records group under, shared by all of them. |
| `SVTYPE` | every record | The type the SV is scored and reported under, read from its first record. On a single-record SV it also names the operation. |
| `OP_TYPE` | every record of a multi-record SV | Named operation from [insilicoSV](https://github.com/PopicLab/insilicoSV). A single-record SV may omit it, or carry the `NA` insilicoSV writes there; either way it goes unused. |
| `SVLEN` | — | When present, the span is `stop = start + abs(SVLEN)`. The absolute value is taken for VCF formats which report negative lengths for deleted segments. Using `SVLEN`, when provided, takes priority over `END` due to pysam autoshifting behavior|
| `TARGET` | dispersed operations | Where the segment is inserted. Absent, an insert lands at its own `stop` — i.e. in tandem. |
| `INSORD` | inserts sharing a `TARGET` | Orders them by insertion order at that position (see below). |
| `TARGET_CHROM` | — | Must equal the record's own chromosome; interchromosomal calls are unsupported. |

Each SV's records are checked before reconstruction:

- every record carries an `SVTYPE`, and every record of a multi-record SV an `OP_TYPE`;
- all records are on one chromosome;
- every `[start, stop)` is non-empty;
- no two intervals overlap (identical spans are fine);
- no `TARGET` lands strictly inside another interval, which would split a segment that interval's
  own operations need; a target exactly on a boundary is fine.

## Procedure

For each SV, svrecon reconstructs the sequence the call implies and asks whether that sequence
exists in the sample.

**1. Reconstruct.** A complex SV may carry several operations, so each SV yields one or more
**queries**: the reference spanning the operations' positions and targets, transformed as the call
describes, plus `--buffer` bp of reference context on each side. If the call is correct, that
query should appear in the sample. A dispersion can produce one or two queries: one if the source
and target buffers overlap and merge into a single window, two independent ones if they don't.
Each SV is reconstructed in isolation, so nearby SVs — called or not — are ignored.

The reconstruction process parses a VCF in [insilicoSV](https://github.com/PopicLab/insilicoSV) format to produce an initial set of disjoint
reference segments, as well as a sequence of deletion, insertion, and inversion operations, where a
complex record can produce more than one. To apply them unambiguously, these operations are ordered:

1. Rightmost-first, so an applied edit never shifts a still-pending operation's coordinates.
2. At a tied position, invert/delete before insert -- they need an untouched segment boundary,
   which a same-position insert would shift.
3. Among tied inserts, by descending INSORD ([insilicoSV](https://github.com/PopicLab/insilicoSV)'s insertion-order field) -- each insert at a
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

- *alignment error* (always) — divergence within the aligned block (mismatches + indels /
  aligned length) ≤ `match_error_threshold`; clips don't count.
- *no large errors* (on by default, `max_contiguous_error: 50`) — minimum contiguous length of a
  failing insertion, deletion, or soft clip.
- *maximum error window* (on by default, `max_window_size: 100` / `max_window_error: 25`) — no
  window of `max_window_size` alignment columns holds `max_window_error` or more error columns
  (mismatches, insertions, deletions, clips)
- *mappable* (on by default, `min_mappable_fraction: 0.95`, `max_unmappable_size: 50`) — inconclusive if less than `min_mappable_fraction` of the query is ACGT, or if it contains a contiguous non-ACGT run of at least `max_unmappable_size`.

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
{
  "svid": "sv16",
  "svtype": "delINV",
  "outcome": "hit",
  "tier": "assembly",
  "query_validations": [
    {
      "query": {
        "chrom": "chr21",
        "svtype": "delINV",
        "grammar": "~AB~->~b~",
        "ref_start": 14952924,
        "ref_end": 14960270,
        "initial_segments": [
          [
            14952924,
            14953024
          ],
          [
            14953024,
            14957080
          ],
          [
            14957080,
            14960170
          ],
          [
            14960170,
            14960270
          ]
        ],
        "ref_segments": [
          [
            14952924,
            14953024
          ],
          [
            14957080,
            14960170
          ],
          [
            14960170,
            14960270
          ]
        ],
        "recon_segments": [
          [
            0,
            100
          ],
          [
            100,
            3190
          ],
          [
            3190,
            3290
          ]
        ],
        "buffer": 100,
        "length": 3290
      },
      "status": "pass",
      "reason": "pass",
      "aligned": true,
      "reference_ambiguous": false,
      "source": "assembly",
      "strand": 1,
      "validations": [
        {
          "source": "assembly",
          "passed": true,
          "status": "pass",
          "reason": "pass",
          "aligned": true,
          "cigar": "3290=",
          "lowest_pass_error": 0.0,
          "lowest_error": 0.0,
          "strand": 1,
          "checks": [
            {
              "status": "pass",
              "detail": "alignment error 0 <= 0.1"
            },
            {
              "status": "pass",
              "detail": "no error run >= 50"
            },
            {
              "status": "pass",
              "detail": "no window >= 25 errors per 100 (worst 0:0)"
            },
            {
              "status": "pass",
              "detail": "mappable 1 >= 0.95, longest nonmappable run 0 < 50"
            }
          ]
        }
      ],
      "ambiguity_validations": []
    }
  ]
}
```

`grammar` is `<reference>-><reconstructed>`: A/B/C... label the SV's distinct intervals in
reference order, `~` is untouched buffer, lowercase means inverted. `~ABC~->~cbC~` (`delINVdup`):
A deleted, C inverted and moved before B, B inverted in place, C also kept in place. One full
grammar string is produced per query -- the `<reference>` half is the same for every query of an
SV, only `<reconstructed>` varies per query.

A multi-operation complex SV produces one `SVValidationResult` with multiple query entries — one
per reconstructed query, not per VCF record. The nesting of the output json mirrors the pipeline: **SV → queries** (one or more reconstructed queries per call, per step 1 above) **→ validations** (one per scorer that ran on that query -- reads, assembly, edlib -- per step 2) **→ checks** (one per configured CIGAR check within that scorer's alignment, per step 3).

## Score summary (`svrecon.summary.json`)

The logged score table also written as
`svrecon.summary.json`, holding the information in two perspectives:
`by_type` is `{svtype: {stat: value}}`, `by_stat` is `{stat: {svtype: value}}`.

## Benchmark

[`workflows/svrecon_benchmark.ipynb`](workflows/svrecon_benchmark.ipynb) measures what the CIGAR
checks buy you. It simulates 1320 SVs per arm (22 types x 20 each x 3 size classes) on hg38 chr21
with [insilicoSV](https://github.com/PopicLab/insilicoSV), then scores two callsets against one assembly, per size class:

| arm | callset | assembly | ideal |
| --- | --- | --- | --- |
| positive | seed 0 | seed 0 | 440/440 hits per size class |
| negative | seed 1 | seed 0 | 0/440 hits per size class |

Each arm is scored under three check configurations, each adding one more check — `similarity`
(bulk alignment error alone; every opt-in check disabled), `similarity--no-large-errors` (adds
`max_contiguous_error: 50`), and `similarity--windowered-errors` (adds the error-window check at
`max_window_size: 100` / `max_window_error: 25`) — so a hit in the negative arm is a false positive
attributable to that configuration:

| arm | mode | small | medium | large | total hits | total inconclusive | precision |
| --- | --- | --- | --- | --- | --- | --- | --- |
| positive | similarity | 440/440 | 440/440 | 439/440 | 1319/1320 | 1/1320 | 1.00 |
| positive | similarity--no-large-errors | 440/440 | 440/440 | 439/440 | 1319/1320 | 1/1320 | 1.00 |
| positive | similarity--windowered-errors | 440/440 | 440/440 | 439/440 | 1319/1320 | 1/1320 | 1.00 |
| negative | similarity | 285/440 | 428/440 | 436/440 | 1149/1320 | 0/1320 | 0.87 |
| negative | similarity--no-large-errors | 43/440 | 0/440 | 0/440 | 43/1320 | 0/1320 | 0.03 |
| negative | similarity--windowered-errors | 0/440 | 0/440 | 0/440 | 0/1320 | 0/1320 | 0.00 |

## Notes

- Outputs go to the `--config` file's directory; without a `--config`, to the working directory.
- Per-chromosome `.mmi` indices are cached (in temp space unless `--chrom-cache` is set) and
  rebuilt automatically when the source FASTA changes.
- Dot plots use wotplot, which accepts only `A`/`C`/`G`/`T`; a query containing `N` is skipped
  unless `--plot-substitute-bases` is given.
