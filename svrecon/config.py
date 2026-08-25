"""CLI config resolution and log-basename derivation."""
import logging
import os
import sys
import tempfile
from pathlib import Path

import yaml

logger = logging.getLogger(__name__)

DEFAULTS = {
    # Main command line defaults
    'config': None,
    'reference': None,
    'sample': None,
    'calls': None,
    'bam': None,
    'classified': None,
    'gap_file': None,
    'n_threads': max(1, os.cpu_count() // 2),
    'verbose': False,

    # Path setting / plotting
    'chrom_cache': None,
    'igv_prefix': '',
    'report': 'none',
    'plot_first_n': 0,
    'plot_aspect': 'auto',
    'plot_substitute_bases': False,
    'plot_title_svid': False,
    'plot_title_location': False,
    'plot_axis_length': False,

    # Shared parameters for validation
    'buffer': 500,
    'auto_buffer_min': 50,        # floor for buffer='auto', also its no-segments fallback
    'auto_buffer_fraction': 0.1,  # buffer='auto' sizes to this fraction of the SV's longest segment
    'match_error_threshold': 0.1,
    'max_indel_size': None,   # if set, fail an alignment carrying an indel at least this long
    'junction_radius': None,  # if set, check the error rate within this radius of every breakpoint

    # Scoring parameters: read-based validation
    'read_error_threshold': 0.1,
    'max_reads_per_site': 1000,

    # Scoring parameters: assembly-based validation (mappy)
    'location_tolerance': float('inf'),
    'assembly_forward_match_only': False,

    # Scoring parameters: edlib fallback validation
    'min_edlib_query': 5000,  # queries at least this long skip the edlib fallback (cost bound)
    # Search-window radius cap for the edlib fallback. This is a COST bound on an
    # unbounded O(n*m) edlib search, kept separate from location_tolerance: the latter
    # only filters mappy hits by genomic position and defaults to infinity (meaningless
    # for unscaffolded, per-contig assemblies whose hit coordinates are contig-local).
    'edlib_fallback_max_tolerance': 10_000_000,

    # Scoring parameters: reference ambiguity check
    'check_reference': False,
}

VALID_PARAM_FNS = {
    # maps valid keys to functions indicating if a given value is valid. A key that
    # doesn't need real value validation (e.g. a path) still needs an entry here --
    # this dict also doubles as the whitelist of settable param names.

    # Main command line defaults
    'config': lambda arg: arg is None or isinstance(arg, str),
    'reference': lambda arg: arg is None or isinstance(arg, str),
    'sample': lambda arg: arg is None or isinstance(arg, str),
    'calls': lambda arg: arg is None or isinstance(arg, str),
    'bam': lambda arg: arg is None or isinstance(arg, str),
    'classified': lambda arg: arg is None or isinstance(arg, str),
    'gap_file': lambda arg: arg is None or isinstance(arg, str),
    'n_threads': lambda arg: isinstance(arg, int),
    'verbose': lambda arg: isinstance(arg, bool),

    # Path setting / plotting
    'chrom_cache': lambda arg: arg is None or isinstance(arg, str),
    'igv_prefix': lambda arg: isinstance(arg, str),
    'report': lambda arg: arg in ('none', 'json'),
    'plot_first_n': lambda arg: isinstance(arg, int),
    'plot_aspect': lambda arg: arg in ('equal', 'auto'),
    'plot_substitute_bases': lambda arg: isinstance(arg, bool),
    'plot_title_svid': lambda arg: isinstance(arg, bool),
    'plot_title_location': lambda arg: isinstance(arg, bool),
    'plot_axis_length': lambda arg: isinstance(arg, bool),

    # Shared parameters for validation
    'buffer': lambda arg: arg == 'auto' or isinstance(arg, int),
    'auto_buffer_min': lambda arg: isinstance(arg, int),
    'auto_buffer_fraction': lambda arg: isinstance(arg, (int, float)),
    'match_error_threshold': lambda arg: isinstance(arg, (int, float)),
    'max_indel_size': lambda arg: arg is None or isinstance(arg, int),
    'junction_radius': lambda arg: arg is None or isinstance(arg, int),

    # Scoring parameters: read-based validation
    'read_error_threshold': lambda arg: isinstance(arg, (int, float)),
    'max_reads_per_site': lambda arg: isinstance(arg, int),

    # Scoring parameters: assembly-based validation (mappy)
    'location_tolerance': lambda arg: isinstance(arg, (int, float)),
    'assembly_forward_match_only': lambda arg: isinstance(arg, bool),

    # Scoring parameters: edlib fallback validation
    'min_edlib_query': lambda arg: isinstance(arg, int),
    'edlib_fallback_max_tolerance': lambda arg: isinstance(arg, int),

    # Scoring parameters: reference ambiguity check
    'check_reference': lambda arg: isinstance(arg, bool),

    # Derived paths (set by Config itself, not user-settable)
    'experiment_dir': lambda arg: isinstance(arg, Path),
    'log_path': lambda arg: isinstance(arg, Path),
    'report_path': lambda arg: arg is None or isinstance(arg, Path),
    'cache_dir': lambda arg: isinstance(arg, Path),
    'img_dir': lambda arg: isinstance(arg, Path),
}

class Config:
    """Merge precedence: CLI flags (`args`) > svrecon --config (YAML) > DEFAULTS."""

    def __init__(self, args):
        self.__dict__.update(DEFAULTS)

        if args.config:
            with open(args.config, "r") as file:
                config_data = yaml.safe_load(file) or {}
            unknown_keys = sorted(set(config_data) - set(VALID_PARAM_FNS))
            if unknown_keys:
                raise ValueError(f"unknown key(s) in {args.config}: {unknown_keys}")
            self.__dict__.update(config_data)

        self.update_from_args(args)

        if self.buffer != 'auto':
            self.buffer = int(self.buffer)

        # Path setup
        self.experiment_dir = Path(args.config).parent.resolve() if args.config else Path.cwd()  # overridden to parent folder of config file if it exists
        self.log_path = self.experiment_dir / f'svrecon.log'
        self.report_path = self.experiment_dir / f'svrecon.report.jsonl' if self.report == 'json' else None
        self.img_dir = self.experiment_dir / 'sv_recon_img'
        if self.plot_first_n:
            self.img_dir.mkdir(parents=True, exist_ok=True)

        # Cache dir for per-chromosome mappy indices
        if self.chrom_cache:
            self.cache_dir = Path(self.chrom_cache)
        else:
            self.cache_dir = Path(tempfile.gettempdir()) / 'mappy_chrom_cache'
        self.cache_dir.mkdir(parents=True, exist_ok=True)

        # Parameter Validation
        for key, value in self.__dict__.items():
            if key not in VALID_PARAM_FNS:
                raise ValueError(f"{key!r} is not a recognized param")
            if not VALID_PARAM_FNS[key](value):
                raise ValueError(f"invalid value for {key!r}: {value!r}")

        # --- logging setup ---
        file_handler = logging.FileHandler(self.log_path, mode='w')
        file_handler.setLevel(logging.INFO)

        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(logging.INFO)

        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s %(levelname)-8s %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S',
            handlers=[file_handler, console_handler],
        )

        logger.info(f'Logging to {self.log_path}')
        if self.config:
            logger.info(f'Loaded svrecon config: {self.config} (outputs -> {self.experiment_dir})')
        logger.info(f'Config:\n{self}')
        logger.info(f"Using {'persistent' if self.chrom_cache else 'temporary'} cache directory: {self.cache_dir}")

    def __repr__(self):
        return '\n'.join(f'  {key}={value!r}' for key, value in vars(self).items())

    def update_from_args(self, args):
        for k, v in vars(args).items():
            if v is None:
                continue
            # Disallow repeat parameters between config and args
            if k in DEFAULTS and getattr(self, k, None) != DEFAULTS[k]:
                raise ValueError(f"'{k}' is set both in the config file ({getattr(self, k)!r}) "
                                 f"and as a CLI flag ({v!r}) -- specify it in only one place")
            setattr(self, k, v)

