import argparse
from collections import defaultdict
import datetime
import hashlib
import logging
import sys
import tempfile
import os

import pandas as pd

from svrecon.align import get_chrom_aligner
from svrecon.config import Config
from svrecon.reads import BamReader
from svrecon.scoring import AlignScorer
from svrecon.util import export_igv_session, load_fasta_to_bytes

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
    parser.add_argument('--groovi_config', help='(Optional) groovi call config used to INFER unset '
                        'params (reference/sample/bam/calls/classified). Was --config previously.',
                        dest='groovi_config', default=None)
    parser.add_argument('--reference', help='Reference genome .fa file', dest='reference', default=None)
    parser.add_argument('--sample', help='Sample genome .fa file(s)', dest='sample', nargs='+', default=None)
    parser.add_argument('--calls', help='VCF containing called SVs', dest='calls', default=None)
    parser.add_argument('--location_tolerance',
                        help='BP tolerance between a mappy hit and the expected SV location; '
                             'default unbounded (recommended for per-contig/unscaffolded assemblies '
                             'whose hit coordinates are contig-local, not genomic)',
                        type=float, default=None)
    parser.add_argument('--buffer', help="Subsequence context buffer, in bp (default 500), or 'auto' to size "
                        "it per-SV to max(100, 10%% of that SV's own longest segment).", dest='buffer', default=None)
    parser.add_argument('--gap_file', help='Tab-delimited file containing regions to omit (e.g., centromere and telomere)',
                        default=None)
    parser.add_argument('--bam', help='BAM file for generating IGV config', dest='bam', default=None)
    parser.add_argument('--classified', help='VCF file of groovi-style classified breakpoints for IGV config',
                        dest='classified', default=None)
    parser.add_argument('--igv_prefix', help='Prefix for igv session paths (default empty)', dest='igv_prefix', default=None)
    parser.add_argument('--chrom_cache', help='Directory to save/load per-chromosome MMI indices.', default=None)
    parser.add_argument('--eval_mode', choices=['assembly', 'reads', 'both', 'none'], default=None,
                        help="Validation source: 'assembly' (default), 'reads' "
                             "(skip the assembly entirely and validate against BAM long reads), "
                             "'both' (assembly first, reads to rescue misses), or 'none' (skip all "
                             "validation -- just reconstruct each SV's alt allele; useful with "
                             "--plot_first_n to get reconstructed-vs-reference dot plots without "
                             "needing a --sample assembly or --bam).")
    parser.add_argument('--read_error_threshold', type=float, default=None,
                        help='Max edlib error rate for a read to validate a reconstruction (read modes; '
                             'default 0.1). Keep <= the assembly error threshold (0.1) for consistent hit/miss calls.')
    parser.add_argument('--min_read_support', type=int, default=None,
                        help='Min number of spanning reads that must clear --read_error_threshold (read modes; default 1).')
    parser.add_argument('--max_reads_per_site', type=int, default=None,
                        help='Cap on candidate reads gathered per SV locus (read modes; default 1000); '
                             'bounds work on deep read pileups.')
    parser.add_argument('--report', choices=['none', 'json'], default=None,
                        help="Write a per-SV evaluation sidecar. 'json' emits "
                             "<logbasename>.eval.jsonl next to the log (one JSON record per SV: "
                             "svid, outcome, tier, and per-segment status/reason/source/error/"
                             "junctions). Default 'none' (log and score table are unchanged).")
    parser.add_argument('--check_reference', action='store_true', default=None,
                        help="Also validate each passing allele against the reference (bulk + junctions); "
                             "if it also validates there (the match isn't specific to the SV -- common "
                             "in repetitive / segmental-dup regions), mark the call inconclusive "
                             "(reason 'reference_match'). Off by default; builds a reference aligner set.")
    parser.add_argument('--junction_window_factor', type=float, default=None,
                        help="Junction-validation window as a multiple of SV size (default 1.5), applied "
                             "to ALL junction checks. Scaling to SV size keeps a real SV's junction signal "
                             "above threshold instead of diluting it in a fixed wide context window.")
    parser.add_argument('--junction_window_min', type=int, default=None,
                        help="Min junction window in bp (default 150). Sets the resolution floor: SVs "
                             "larger than ~2*min*T/(1-T) are resolvable against the reference.")
    parser.add_argument('--junction_window_max', type=int, default=None,
                        help="Max junction window in bp (default 300); caps the window for large SVs.")
    parser.add_argument('--plot_first_n', type=int, default=None,
                        help="Write a dot-plot PNG (reconstructed subsequence vs. the validating real-data "
                             "sequence) for the first N SV calls of each SV type, into <output_dir>/img/ "
                             "(default 0 -- no plots).")
    parser.add_argument('--plot_aspect', choices=['equal', 'auto'], default=None,
                        help="Dot-plot axes aspect (default 'auto'): 'auto' keeps the plot square; "
                             "'equal' is true-to-scale (1bp=1bp) but can squeeze asymmetric SVs into a sliver.")
    parser.add_argument('--assembly_forward_match_only', action='store_true', default=None,
                        help="Only accept forward-strand alignments during assembly validation "
                             "(reverse-strand hits are filtered out before the error/junction checks). "
                             "Off by default (maps to both strands).")
    args = parser.parse_args()

    config = Config(args)

    eval_mode = config.eval_mode

    logger.info('Initializing scorer and loading callset')
    scorer = AlignScorer(config.calls, config.buffer, config.gap_file)
    scorer.eval_mode = eval_mode
    scorer.read_error_threshold = config.read_error_threshold
    scorer.min_read_support = config.min_read_support
    scorer.max_reads_per_site = config.max_reads_per_site
    scorer.check_reference = config.check_reference
    scorer.junction_window_factor = config.junction_window_factor
    scorer.junction_window_min = config.junction_window_min
    scorer.junction_window_max = config.junction_window_max
    scorer.assembly_forward_match_only = config.assembly_forward_match_only

    logger.info('Finding relevant chromosomes')
    chroms = set()
    for records in scorer.variants.values():
        for record in records:
            chroms.add(record.chrom)
            if 'TARGET_CHROM' in record.info:
                chroms.add(record.info['TARGET_CHROM'])
    logger.info(f'Found {len(chroms)} referenced chromosomes in callset')

    logger.info('Loading reference bytearrays')
    scorer.ref = load_fasta_to_bytes(config.reference, chroms)
    logger.info(f'Loaded {len(scorer.ref)} reference chromosomes.')

    scorer.sample = []
    scorer.aligners = defaultdict(list)
    scorer.reference_aligners = None

    # Aligner build params + cache dir, shared by the assembly aligners and (when
    # --check_reference is on) the reference aligners.
    align_params = {
        'preset': 'map-hifi',
        'k': 15,
        'w': 5,
        'best_n': 100,
        'min_cnt': 1,
        'min_dp_score': 10,
        'min_chain_score': 1,
    }
    if config.chrom_cache:
        cache_dir = config.chrom_cache
        logger.info(f'Using persistent cache directory: {cache_dir}')
    else:
        cache_dir = os.path.join(tempfile.gettempdir(), 'mappy_chrom_cache')
        logger.info(f'Using temporary cache directory: {cache_dir}')
    os.makedirs(cache_dir, exist_ok=True)

    def build_chrom_aligners(fasta, into, label):
        """Build/load per-chromosome mappy aligners for `fasta` into the `into` dict."""
        path_hash = hashlib.md5(os.path.abspath(fasta).encode('utf-8')).hexdigest()[:8]
        fa_cache = os.path.join(cache_dir, f'{os.path.basename(fasta)}_{path_hash}')
        os.makedirs(fa_cache, exist_ok=True)
        for chrom in chroms:
            aligner = get_chrom_aligner(fasta, chrom, fa_cache, align_params, threads=32)
            if aligner:
                into[chrom].append(aligner)
            else:
                logger.warning(f'No sequences found for {chrom} in {label} FASTA {fasta}.')

    if eval_mode in ('assembly', 'both'):
        logger.info('Loading sample assembly bytearrays')
        scorer.sample = [load_fasta_to_bytes(samp, chroms) for samp in config.sample]
        logger.info(f'Initialized {len(scorer.sample)} sample assembly file(s).')
        logger.info('Building/loading per-chromosome aligners...')
        for samp in config.sample:
            build_chrom_aligners(samp, scorer.aligners, 'sample')
        logger.info('Aligners ready')
    else:
        logger.info(f"eval_mode={eval_mode!r}: skipping sample assembly load and aligner build.")

    if config.check_reference:
        logger.info('Building/loading reference aligner(s) for --check_reference (mappy)...')
        scorer.reference_aligners = defaultdict(list)
        build_chrom_aligners(config.reference, scorer.reference_aligners, 'reference')
        logger.info('Reference aligners ready')

    if eval_mode in ('reads', 'both'):
        if not config.bam:
            logger.error('--eval_mode reads/both requires a BAM (--bam) for read-based evaluation.')
            sys.exit(1)
        logger.info(f'Initializing thread-safe BAM reader: {config.bam}')
        scorer.bam_reader = BamReader(config.bam)

    img_dir = config.experiment_dir / 'sv_recon_img'

    precision, correct_calls, total_calls, inconclusive_calls, skipped_calls, assembly_hits, read_hits = scorer.score_all(
        location_tolerance=config.location_tolerance, report_path=config.report_path,
        plot_first_n=config.plot_first_n, plot_out_dir=img_dir, plot_aspect=config.plot_aspect)

    df = pd.DataFrame({
        'correct_calls': correct_calls,
        'total_calls': total_calls,
        'inconclusive': inconclusive_calls,
        'skipped': skipped_calls,
        'precision': precision,
        'assembly_hits': assembly_hits,
        'read_hits': read_hits,
    })
    # count columns: missing SV-type keys (e.g. a type with 0 hits) -> 0, not NaN
    for col in ('correct_calls', 'total_calls', 'inconclusive', 'skipped', 'assembly_hits', 'read_hits'):
        df[col] = df[col].fillna(0).astype(int)
    df.sort_values('total_calls', ascending=False, inplace=True)

    # to_string() prints all columns (default repr truncates the middle ones)
    logger.info('Score table:\n' + df.to_string())

    export_igv_session(config.calls, config.bam, config.classified, timestamp, config.igv_prefix)


if __name__ == '__main__':
    main()
