import argparse
import datetime
import json
import logging

import pandas as pd

from svrecon.config import Config
from svrecon.scoring import CallsetScorer
from svrecon.utils import export_igv_session

logger = logging.getLogger(__name__)


def main():
    timestamp = datetime.datetime.now().strftime('%Y-%m-%d_%H%M%S')

    parser = argparse.ArgumentParser(prog='svrecon',
                                     description='Score VCF SV calls against a reference and sample genome')
    # An svrecon config (YAML whose keys mirror these flags) supplies any of the params below;
    # logs and reports are written to that config's directory. Explicit CLI flags override the
    # config. Mergeable params default to None so we can tell "unset" from an explicit value.
    parser.add_argument('--config', help='svrecon YAML config; keys mirror these flags. Logs and '
                        'reports are written to this file\'s directory.', dest='config')
    parser.add_argument('--reference', help='Reference genome .fa file', dest='reference', default=None)
    parser.add_argument('--sample', help='Sample genome .fa file(s). Presence enables assembly-based '
                        'validation.', dest='sample', default=None)
    parser.add_argument('--calls', help='VCF containing called SVs', dest='calls', default=None)
    parser.add_argument('--location-tolerance',
                        help='BP tolerance between a mappy hit and the expected SV location; '
                             'default unbounded (recommended for per-contig/unscaffolded assemblies '
                             'whose hit coordinates are contig-local, not genomic)',
                        type=float, default=None)
    parser.add_argument('--buffer', help="Query context buffer, in bp (default 500), or 'auto' to size "
                        "it per-SV to max(50, 10%% of that SV's own longest segment).", dest='buffer', default=None)
    parser.add_argument('--gap-file', help='Tab-delimited file containing regions to omit (e.g., centromere and telomere)',
                        default=None)
    parser.add_argument('--bam', help='BAM file. Presence enables read-based validation; also used '
                        'for generating IGV config.', dest='bam', default=None)
    parser.add_argument('--classified', help='VCF file of groovi-style classified breakpoints for IGV config',
                        dest='classified', default=None)
    parser.add_argument('--igv-prefix', help='Prefix for igv session paths (default empty)', dest='igv_prefix', default=None)
    parser.add_argument('--chrom-cache', help='Directory to save/load per-chromosome MMI indices.', default=None)
    parser.add_argument('--read-error-threshold', type=float, default=None,
                        help='Max edlib error rate for a read to validate a reconstruction (read modes; '
                             'default 0.1). Keep <= the assembly error threshold (0.1) for consistent hit/miss calls.')
    parser.add_argument('--match-error-threshold', type=float, default=None,
                        help='Max error rate for an assembly/edlib/reference alignment to validate a '
                             'reconstruction (default 0.1).')
    parser.add_argument('--n-threads', type=int, default=None,
                        help='Worker threads for scoring SVs (default: half the CPU count).')
    parser.add_argument('--verbose', action='store_true', default=None,
                        help='Log more detail per SV.')
    parser.add_argument('--max-reads-per-site', type=int, default=None,
                        help='Cap on candidate reads gathered per SV locus (read modes; default 1000); '
                             'bounds work on deep read pileups.')
    parser.add_argument('--report', choices=['none', 'json'], default=None,
                        help="Write a per-SV evaluation sidecar. 'json' emits svrecon.report.jsonl "
                             "next to the log (one JSON record per SV: svid, outcome, tier, and "
                             "each query's status/reason/source/errors/checks). Default 'none' "
                             "(log and score table are unchanged).")
    parser.add_argument('--check-reference', action='store_true', default=None,
                        help="Also validate each passing allele against the reference; if it validates "
                             "there too (the alignment isn't specific to the SV -- common in repetitive "
                             "/ segmental-dup regions), mark the call inconclusive "
                             "(reference_ambiguous). Off by default; builds a reference aligner set.")
    parser.add_argument('--plot-first-n', type=int, default=None,
                        help="Write a dot-plot PNG (reconstructed subsequence vs. the validating real-data "
                             "sequence) for the first N SV calls of each SV type, into <output_dir>/img/ "
                             "(default 0 -- no plots).")
    parser.add_argument('--plot-substitute-bases', action='store_true', default=None,
                        help='Replace non-ACGT bases with random ACGT bases when plotting, instead '
                             'of skipping the plot. Affects plots only, never validation.')
    parser.add_argument('--plot-aspect', choices=['equal', 'auto'], default=None,
                        help="Dot-plot axes aspect (default 'auto'): 'auto' keeps the plot square; "
                             "'equal' is true-to-scale (1bp=1bp) but can squeeze asymmetric SVs into a sliver.")
    parser.add_argument('--plot-title-svid', action='store_true', default=None,
                        help='Include the SV id in the plot title.')
    parser.add_argument('--plot-title-location', action='store_true', default=None,
                        help='Include the locus in the plot title.')
    parser.add_argument('--plot-axis-length', action='store_true', default=None,
                        help='Include each sequence length in its axis label.')
    parser.add_argument('--assembly-forward-match-only', action='store_true', default=None,
                        help="Only accept forward-strand alignments during assembly validation "
                             "(reverse-strand hits are filtered out before the error/segment checks). "
                             "Off by default (maps to both strands).")
    args = parser.parse_args()

    logger.info('Initializing scorer and loading callset')
    config = Config(args)
    scorer = CallsetScorer(config)

    precision, correct_calls, total_calls, inconclusive_calls, assembly_hits, read_hits = scorer.score_all()

    df = pd.DataFrame({
        'correct_calls': correct_calls,
        'total_calls': total_calls,
        'inconclusive': inconclusive_calls,
        'precision': precision,
        'assembly_hits': assembly_hits,
        'read_hits': read_hits,
    })
    # count columns: missing SV-type keys (e.g. a type with 0 hits) -> 0, not NaN
    for col in ('correct_calls', 'total_calls', 'inconclusive', 'assembly_hits', 'read_hits'):
        df[col] = df[col].fillna(0).astype(int)
    df.sort_values('total_calls', ascending=False, inplace=True)

    # to_string() prints all columns (default repr truncates the middle ones)
    logger.info('Score table:\n' + df.to_string())

    summary = {'by_stat': df.to_dict(), 'by_type': df.to_dict(orient='index')}
    with open(config.summary_json_path, 'w') as f:
        json.dump(summary, f, indent=2)
    logger.info(f'Wrote per-type score summary: {config.summary_json_path}')

    export_igv_session(config.calls, config.bam, config.classified, timestamp, config.igv_prefix)


if __name__ == '__main__':
    main()
