"""svrecon command-line entry point."""
import argparse
import datetime
import hashlib
import logging
import os
import sys
import tempfile
from collections import defaultdict

import pandas as pd
import yaml

from svrecon.align import get_chrom_aligner
from svrecon.config import log_name_from_config, update_args_from_groovi_config
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
    parser.add_argument('--buffer', help='Subsequence context buffer (default 500)', type=int, dest='buffer', default=None)
    parser.add_argument('--gap_file', help='Tab-delimited file containing regions to omit (e.g., centromere and telomere)',
                        default=None)
    parser.add_argument('--bam', help='BAM file for generating IGV config', dest='bam', default=None)
    parser.add_argument('--classified', help='VCF file of groovi-style classified breakpoints for IGV config',
                        dest='classified', default=None)
    parser.add_argument('--igv_prefix', help='Prefix for igv session paths (default empty)', dest='igv_prefix', default=None)
    parser.add_argument('--chrom_cache', help='Directory to save/load per-chromosome MMI indices.', default=None)
    parser.add_argument('--eval_mode', choices=['assembly', 'reads', 'both'], default=None,
                        help="Validation source: 'assembly' (default), 'reads' "
                             "(skip the assembly entirely and validate against BAM long reads), or "
                             "'both' (assembly first, reads to rescue misses).")
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
    args = parser.parse_args()

    # --- Resolve the effective config: CLI flag > svrecon --config value > groovi inference > default ---
    MERGE_KEYS = ['reference', 'sample', 'calls', 'bam', 'classified', 'gap_file', 'chrom_cache',
                  'igv_prefix', 'eval_mode', 'buffer', 'location_tolerance', 'read_error_threshold',
                  'min_read_support', 'max_reads_per_site', 'report', 'check_reference',
                  'junction_window_factor', 'junction_window_min',
                  'junction_window_max', 'groovi_config']
    cfg, unknown_keys = {}, []
    if args.config:
        with open(args.config) as f:
            cfg = yaml.safe_load(f) or {}
        unknown_keys = sorted(set(cfg) - set(MERGE_KEYS))
        for k in MERGE_KEYS:                       # svrecon config fills anything not set on the CLI
            if getattr(args, k) is None and k in cfg:
                setattr(args, k, cfg[k])

    if args.groovi_config:                          # groovi inference fills anything still unset
        args = update_args_from_groovi_config(args, args.groovi_config)

    defaults = {'eval_mode': 'assembly', 'buffer': 500, 'location_tolerance': float('inf'),
                'igv_prefix': '', 'read_error_threshold': 0.1, 'min_read_support': 1,
                'max_reads_per_site': 1000, 'report': 'none', 'check_reference': False,
                'junction_window_factor': 1.5, 'junction_window_min': 150,
                'junction_window_max': 300}
    for k, v in defaults.items():
        if getattr(args, k) is None:
            setattr(args, k, v)

    # YAML already infers int/float for numeric params (write infinity as `.inf`); a bad value
    # crashes on use. But a wrong eval_mode/report from the config does NOT crash -- it silently
    # mis-scores or skips the report -- so validate those. Also accept a scalar `sample:` path.
    if isinstance(args.sample, str):
        args.sample = [args.sample]
    if args.eval_mode not in ('assembly', 'reads', 'both'):
        parser.error(f"eval_mode must be assembly|reads|both, got {args.eval_mode!r}")
    if args.report not in ('none', 'json'):
        parser.error(f"report must be none|json, got {args.report!r}")

    # Outputs live beside the svrecon --config: the experiment directory identifies the run, so
    # fixed names (a re-run in the same dir overwrites -- use separate configs/dirs to compare
    # modes). Without a --config, fall back to ./logs with a run-specific name (the groovi
    # experiment name, else a timestamp) so ad-hoc runs don't clobber each other.
    if args.config:
        output_dir = os.path.dirname(os.path.abspath(args.config))
        base = 'svrecon'
    else:
        output_dir = './logs'
        base = log_name_from_config(args.groovi_config) or f'svrecon_{timestamp}'
    os.makedirs(output_dir, exist_ok=True)
    log_filename = os.path.join(output_dir, f'{base}.log')
    report_path = os.path.join(output_dir, f'{base}.report.jsonl') if args.report == 'json' else None

    file_handler = logging.FileHandler(log_filename, mode='w')
    file_handler.setLevel(logging.DEBUG)

    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)

    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s %(levelname)-8s %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        handlers=[file_handler, console_handler]
    )

    logger.info(f'Logging to {log_filename}')
    if args.config:
        logger.info(f'Loaded svrecon config: {args.config} (outputs -> {output_dir})')
    if unknown_keys:
        logger.warning(f'Ignoring unrecognized keys in {args.config}: {unknown_keys}')
    logger.info(f'Config: {vars(args)}')


    eval_mode = args.eval_mode

    logger.info('Initializing scorer and loading callset')
    scorer = AlignScorer(args.calls, args.buffer, args.gap_file)
    scorer.eval_mode = eval_mode
    scorer.read_error_threshold = args.read_error_threshold
    scorer.min_read_support = args.min_read_support
    scorer.max_reads_per_site = args.max_reads_per_site
    scorer.check_reference = args.check_reference
    scorer.junction_window_factor = args.junction_window_factor
    scorer.junction_window_min = args.junction_window_min
    scorer.junction_window_max = args.junction_window_max

    logger.info('Finding relevant chromosomes')
    chroms = set()
    for records in scorer.variants.values():
        for record in records:
            chroms.add(record.chrom)
            if 'TARGET_CHROM' in record.info:
                chroms.add(record.info['TARGET_CHROM'])
    logger.info(f'Found {len(chroms)} referenced chromosomes in callset')

    logger.info('Loading reference bytearrays')
    scorer.ref = load_fasta_to_bytes(args.reference, chroms)
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
    if args.chrom_cache:
        cache_dir = args.chrom_cache
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
        scorer.sample = [load_fasta_to_bytes(samp, chroms) for samp in args.sample]
        logger.info(f'Initialized {len(scorer.sample)} sample assembly file(s).')
        logger.info('Building/loading per-chromosome aligners...')
        for samp in args.sample:
            build_chrom_aligners(samp, scorer.aligners, 'sample')
        logger.info('Aligners ready')
    else:
        logger.info('Read-only eval mode: skipping sample assembly load and aligner build.')

    if args.check_reference:
        logger.info('Building/loading reference aligner(s) for --check_reference (mappy)...')
        scorer.reference_aligners = defaultdict(list)
        build_chrom_aligners(args.reference, scorer.reference_aligners, 'reference')
        logger.info('Reference aligners ready')

    if eval_mode in ('reads', 'both'):
        if not args.bam:
            logger.error('--eval_mode reads/both requires a BAM (--bam) for read-based evaluation.')
            sys.exit(1)
        logger.info(f'Initializing thread-safe BAM reader: {args.bam}')
        scorer.bam_reader = BamReader(args.bam)

    precision, correct_calls, total_calls, inconclusive_calls, assembly_hits, read_hits = scorer.score_all(
        location_tolerance=args.location_tolerance, report_path=report_path)

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

    export_igv_session(args.calls, args.bam, args.classified, timestamp, args.igv_prefix)


if __name__ == '__main__':
    main()
