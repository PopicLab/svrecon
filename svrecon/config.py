"""CLI config resolution: groovi-config inference and log-basename derivation."""
import logging
import os
from pathlib import Path

import yaml

logger = logging.getLogger(__name__)


def update_args_from_groovi_config(args, groovi_config_path):
    """Infer parameters from a groovi call config file.

    Fills only params that are still unset (None) on `args`, so it never overrides a
    value given explicitly on the CLI or in an svrecon --config. Path resolution is
    sensitive to the internal groovi folder layout (see README notes)."""
    with open(groovi_config_path, 'r') as file:
        config_data = yaml.safe_load(file)

        experiment_dir = str(Path(groovi_config_path).parent.resolve())
        results_dir = os.path.join(experiment_dir, "results")
        if args.calls is None:
            args.calls = os.path.join(results_dir, 'groovi.vcf')

        exp_path = Path(experiment_dir)
        prefix = next(p for p in exp_path.parents if p.name == 'data').parent
        bam_path = Path(config_data['bam'])
        if not bam_path:
            data_index = bam_path.parts.index('data')
            bam_path = Path(*bam_path.parts[data_index:])

        bam = prefix / bam_path

        if not args.bam:
            args.bam = str(bam)

        sample = bam.parent.parent / 'VCF/sim.fa'

        if sample.exists() and args.sample is None:
            args.sample = str(sample)

        fa_path = Path(config_data['fa'])
        data_index = fa_path.parts.index('data')
        fa_path = Path(*fa_path.parts[data_index:])

        if args.classified is None:
            args.classified = os.path.join(experiment_dir, "results/groovi_bkps_classified.vcf")

        if args.reference is None:
            args.reference = str(prefix / fa_path)

        logger.info(f'Config updated: {vars(args)}')
        logger.debug(f'Full groovi config: {config_data}')
    return args


def log_name_from_config(config):
    """Derive a stable log basename from the config's experiment directory.

    The experiment is identified by the config's directory relative to an
    `experiments/` ancestor, with path components joined by dots. For example:
        .../experiments/HG00733/groovi.yaml            -> svrecon.HG00733
        .../experiments/hg002/latest_filtering/groovi.yaml
                                                       -> svrecon.hg002.latest_filtering
    Falls back to the config's parent directory name when no `experiments/`
    ancestor is present, and returns None if no config is given.
    """
    if not config:
        return None

    exp_dir = os.path.dirname(os.path.abspath(config))
    parts = exp_dir.split(os.sep)
    if 'experiments' in parts:
        idx = parts.index('experiments')
        rel_parts = parts[idx + 1:]
    else:
        rel_parts = parts[-1:]

    rel_parts = [p for p in rel_parts if p]
    if not rel_parts:
        return None
    return 'svrecon.' + '.'.join(rel_parts)
