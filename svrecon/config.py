"""CLI config resolution: groovi-config inference and log-basename derivation."""
import datetime
import logging
import sys
from pathlib import Path

import yaml

logger = logging.getLogger(__name__)

DEFAULTS = {
    'config': None,
    'reference': None,
    'sample': None,
    'calls': None,
    'bam': None,
    'classified': None,
    'gap_file': None,
    'chrom_cache': None,
    'groovi_config': None,
    'eval_mode': 'assembly',
    'buffer': 500,
    'location_tolerance': float('inf'),
    'igv_prefix': '',
    'read_error_threshold': 0.1,
    'min_read_support': 1,
    'max_reads_per_site': 1000,
    'report': 'none',
    'check_reference': False,
    'junction_window_factor': 1.5,
    'junction_window_min': 150,
    'junction_window_max': 300,
    'plot_first_n': 0,
    'plot_aspect': 'auto',
    'assembly_forward_match_only': False,
}

VALID_PARAM_FNS = {
    # maps valid keys to functions indicating if a given value is valid. A key that
    # doesn't need real value validation (e.g. a path) still needs an entry here --
    # this dict also doubles as the whitelist of settable param names.
    'config': lambda arg: arg is None or isinstance(arg, str),
    'reference': lambda arg: arg is None or isinstance(arg, str),
    'sample': lambda arg: arg is None or isinstance(arg, (str, list)),
    'calls': lambda arg: arg is None or isinstance(arg, str),
    'bam': lambda arg: arg is None or isinstance(arg, str),
    'classified': lambda arg: arg is None or isinstance(arg, str),
    'gap_file': lambda arg: arg is None or isinstance(arg, str),
    'chrom_cache': lambda arg: arg is None or isinstance(arg, str),
    'groovi_config': lambda arg: arg is None or isinstance(arg, str),
    'igv_prefix': lambda arg: isinstance(arg, str),
    'eval_mode': lambda arg: arg in ('assembly', 'reads', 'both', 'none'),
    'report': lambda arg: arg in ('none', 'json'),
    'plot_aspect': lambda arg: arg in ('equal', 'auto'),
    'buffer': lambda arg: arg == 'auto' or isinstance(arg, int),
    'location_tolerance': lambda arg: isinstance(arg, (int, float)),
    'read_error_threshold': lambda arg: isinstance(arg, (int, float)),
    'min_read_support': lambda arg: isinstance(arg, int),
    'max_reads_per_site': lambda arg: isinstance(arg, int),
    'junction_window_factor': lambda arg: isinstance(arg, (int, float)),
    'junction_window_min': lambda arg: isinstance(arg, int),
    'junction_window_max': lambda arg: isinstance(arg, int),
    'plot_first_n': lambda arg: isinstance(arg, int),
    'check_reference': lambda arg: isinstance(arg, bool),
    'assembly_forward_match_only': lambda arg: isinstance(arg, bool),
    'experiment_dir': lambda arg: isinstance(arg, Path),
    'log_path': lambda arg: isinstance(arg, Path),
    'report_path': lambda arg: arg is None or isinstance(arg, Path),
}

class Config:
    """Merge precedence: CLI flags (`args`) > svrecon --config (YAML) > groovi inference > DEFAULTS."""

    def __init__(self, args):
        timestamp = datetime.datetime.now().strftime('%Y-%m-%d_%H%M%S')
        logging.info(f"initiating svrecon at {timestamp}")

        self.__dict__.update(DEFAULTS)
        self.experiment_dir = Path.cwd()  # overridden to parent folder of config file if it exists

        if args.config:
            with open(args.config, "r") as file:
                config_data = yaml.safe_load(file) or {}
            unknown_keys = sorted(set(config_data) - set(VALID_PARAM_FNS))
            if unknown_keys:
                raise ValueError(f"unknown key(s) in {args.config}: {unknown_keys}")
            self.__dict__.update(config_data)
            self.experiment_dir = Path(args.config).parent.resolve()

        self.update_from_args(args)

        if isinstance(self.sample, str):
            self.sample = [self.sample]
        if self.buffer != 'auto':
            self.buffer = int(self.buffer)

        base = 'svrecon'
        self.log_path = self.experiment_dir / f'{base}.log'
        self.report_path = self.experiment_dir / f'{base}.report.jsonl' if self.report == 'json' else None

        # Parameter Validation
        for key, value in self.__dict__.items():
            if key not in VALID_PARAM_FNS:
                raise ValueError(f"{key!r} is not a recognized param")
            if not VALID_PARAM_FNS[key](value):
                raise ValueError(f"invalid value for {key!r}: {value!r}")

        # --- logging setup ---
        file_handler = logging.FileHandler(self.log_path, mode='w')
        file_handler.setLevel(logging.DEBUG)

        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setLevel(logging.INFO)

        logging.basicConfig(
            level=logging.DEBUG,
            format='%(asctime)s %(levelname)-8s %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S',
            handlers=[file_handler, console_handler],
        )

        logger.info(f'Logging to {self.log_path}')
        if self.config:
            logger.info(f'Loaded svrecon config: {self.config} (outputs -> {self.experiment_dir})')
        logger.info(f'Config: {vars(self)}')

    def update_from_args(self, args):
        for k, v in vars(args).items():
            if v is None:
                continue
            # Disallow repeat parameters between config and args
            if k in DEFAULTS and getattr(self, k, None) != DEFAULTS[k]:
                raise ValueError(f"'{k}' is set both in the config file ({getattr(self, k)!r}) "
                                 f"and as a CLI flag ({v!r}) -- specify it in only one place")
            setattr(self, k, v)

