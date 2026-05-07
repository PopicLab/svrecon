import logging
import os
from pathlib import Path
from typing import Dict, List, Tuple, Union

import edlib
import mappy
import pysam
import yaml
from pysam import VariantRecord

logger = logging.getLogger(__name__)


def get_chrom_aligner(sample_fasta: str, chrom: str, cache_dir: str, align_params: dict,
                      threads: int = 4) -> mappy.Aligner:
    """Loads a per-chromosome MMI index, building it first if it doesn't exist."""
    mmi_path = os.path.join(cache_dir, f"{chrom}.mmi")

    if os.path.exists(mmi_path):
        logger.info(f"Loading existing index for {chrom} from {mmi_path}")
        return mappy.Aligner(mmi_path, **align_params)

    logger.info(f"Index not found for {chrom}. Extracting sequences...")
    chrom_fa_path = os.path.join(cache_dir, f"{chrom}.fa")
    written_contigs = 0

    with pysam.FastaFile(sample_fasta) as fasta:
        with open(chrom_fa_path, 'w') as out_f:
            for ref in fasta.references:
                if ref == chrom or ref.startswith(f"{chrom}_"):
                    seq = fasta.fetch(ref)
                    out_f.write(f">{ref}\n{seq}\n")
                    written_contigs += 1

    if written_contigs == 0:
        if os.path.exists(chrom_fa_path):
            os.remove(chrom_fa_path)
        return None

    logger.info(f"Building MMI index for {chrom} at {mmi_path}...")
    aligner = mappy.Aligner(chrom_fa_path, n_threads=threads, fn_idx_out=mmi_path, **align_params)

    if os.path.exists(chrom_fa_path):
        os.remove(chrom_fa_path)

    return aligner


def update_args_from_config(args):
    """Infer parameters from groovi config files"""
    with open(args.config, 'r') as file:
        config_data = yaml.safe_load(file)

        experiment_dir = str(Path(args.config).parent.resolve())
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

        args.classified = os.path.join(experiment_dir, "results/groovi_bkps_classified.vcf")

        if args.reference is None:
            args.reference = str(prefix / fa_path)

        logger.info(f'Config updated: {vars(args)}')
    return args


def get_start_stop(rec: VariantRecord) -> Tuple[int, int]:
    """Converts VCF start/stop coordinates to Python style"""
    start = rec.start - 1
    stop = rec.stop
    return start, stop


def reverse_complement(seq: Union[str, List[str]]) -> List[str]:
    """Returns the reverse complement of a given DNA sequence string."""
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
                  'a': 't', 'c': 'g', 'g': 'c', 't': 'a',
                  'N': 'N', 'n': 'n'}
    return [complement[x] for x in reversed(seq)]


def check_match(alignment, query):
    """Calculate the independent error rate of a mappy alignment with the query sequence"""
    aligned_query_segment_length = alignment.q_en - alignment.q_st
    len_unaligned = len(query) - aligned_query_segment_length

    nm_total = alignment.NM + len_unaligned
    len_norm = alignment.blen + len_unaligned

    return nm_total / len_norm if len_norm > 0 else 1.0


def load_fasta_to_bytes(filename: str, chroms) -> Dict[str, bytearray]:
    """Loads target sequences into a dictionary of bytearrays."""
    data = {}
    with pysam.FastaFile(filename) as f:
        fasta_refs = f.references
        for chrom in chroms:
            chrom_lower = chrom.lower()
            for ref in fasta_refs:
                ref_lower = ref.lower()
                if ref_lower == chrom_lower or ref_lower.startswith(f"{chrom_lower}_"):
                    sequence_string = f.fetch(reference=ref)
                    data[ref] = bytearray(sequence_string, 'ascii')
    return data


def run_edlib_fallback(query_seq: str, chrom: str, location: int, initial_buffer: int, max_tolerance: int,
                       error_threshold: float, sample_dict: Dict[str, bytearray]) -> float:
    """
    Executes an exact edit distance search via edlib, exponentially expanding the search window
    starting from initial_buffer up to max_tolerance.
    Returns the calculated error rate, exiting early if it falls below the error_threshold.
    """
    if sample_dict is None:
        return 1.0

    target_keys = [k for k in sample_dict.keys() if k == chrom or k.startswith(f"{chrom}_")]
    if not target_keys:
        return 1.0

    best_error = 1.0
    current_tolerance = initial_buffer

    prev_bounds = {k: (-1, -1) for k in target_keys}

    while True:
        expanded_any = False

        for target_key in target_keys:
            target_len = len(sample_dict[target_key])

            search_start = max(0, location - current_tolerance)
            search_end = min(target_len, location + len(query_seq) + current_tolerance)

            if (search_start, search_end) == prev_bounds[target_key]:
                continue

            expanded_any = True
            prev_bounds[target_key] = (search_start, search_end)

            target_seq = sample_dict[target_key][search_start:search_end].decode('ascii')

            result = edlib.align(query_seq.upper(), target_seq.upper(), mode="HW", task="path")

            if result and result['editDistance'] >= 0:
                edit_dist = result['editDistance']

                if result.get('locations') and result['locations'][0][0] is not None and result['locations'][0][1] is not None:
                    loc = result['locations'][0]
                    target_match_len = loc[1] - loc[0] + 1
                    denominator = max(len(query_seq), target_match_len)
                else:
                    denominator = len(query_seq)

                error_rate = edit_dist / denominator if denominator > 0 else 1.0
                best_error = min(best_error, error_rate)

                if best_error <= error_threshold:
                    return best_error

        if not expanded_any or current_tolerance >= max_tolerance:
            break

        current_tolerance = max(1, current_tolerance * 10)
        if current_tolerance > max_tolerance:
            current_tolerance = max_tolerance

    return best_error