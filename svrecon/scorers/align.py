"""Alignment + error-rate primitives: mappy aligner caching, segment validation, edlib scoring."""
import hashlib
import logging
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Tuple, Union

import edlib
import mappy
import pysam

from svrecon.scorers.utils import Cigar

logger = logging.getLogger(__name__)


@dataclass
class SeqSegmentValidationResult:
    error: float
    passed: bool

    def jsonify(self) -> dict:
        return {'error': self.error, 'passed': self.passed}


def validate_segments_from_cigar(cigar: Cigar, segments: List[Tuple[int, int]],
                                 error_threshold: float = 0.1) -> List[SeqSegmentValidationResult]:
    """Scores each segment [start, end) against the same error threshold; one result per
    segment, verdict = all(r.passed). Segments tile the query (reconstruct.py), and at a
    shared boundary i of [a, i) and [i, b), a deletion anchored at i counts toward
    [i, b) -- so a deletion is checked implicitly through its flanking segments.
    Clipped bases are error, so an unaligned segment fails at 1.0."""
    rates = [cigar.get_window_error_rate(start, end) for start, end in segments]
    return [SeqSegmentValidationResult(error=float(f'{rate:.4g}'), passed=rate <= error_threshold)
            for rate in rates]


def get_chrom_aligner(sample_fasta: str, chrom: str, cache_dir: str, align_params: dict,
                      threads: int = 4) -> mappy.Aligner:
    """Loads a per-chromosome MMI index, building it first if it doesn't exist.
    Raises ValueError if the FASTA has no sequences for the chromosome."""
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
        raise ValueError(f"no sequences found for chromosome {chrom!r} in {sample_fasta}")

    logger.info(f"Building MMI index for {chrom} at {mmi_path}...")
    aligner = mappy.Aligner(chrom_fa_path, n_threads=threads, fn_idx_out=mmi_path, **align_params)

    if os.path.exists(chrom_fa_path):
        os.remove(chrom_fa_path)

    return aligner


# Aligner build params, shared by the assembly aligners and (when
# --check_reference is on) the reference aligners.
ALIGN_PARAMS = {
    'preset': 'map-hifi',
    'k': 15,
    'w': 5,
    'best_n': 100,
    'min_cnt': 1,
    'min_dp_score': 10,
    'min_chain_score': 1,
}


def build_chrom_aligners(fasta_path: str, chroms: Iterable[str], cache_dir: Path) -> Dict[str, mappy.Aligner]:
    """Build/load per-chromosome mappy aligners for `fasta_path`, keyed by chromosome."""
    path_hash = hashlib.md5(str(Path(fasta_path).resolve()).encode('utf-8')).hexdigest()[:8]
    fa_cache = cache_dir / f'{Path(fasta_path).name}_{path_hash}'
    fa_cache.mkdir(parents=True, exist_ok=True)
    return {chrom: get_chrom_aligner(fasta_path, chrom, str(fa_cache), ALIGN_PARAMS, threads=32)
            for chrom in chroms}


# TODO: rename to something like get normalized score, change query param from str to query object
def check_match(alignment, query):
    """Calculate the independent error rate of a mappy alignment with the query sequence"""
    aligned_query_segment_length = alignment.q_en - alignment.q_st
    len_unaligned = len(query) - aligned_query_segment_length

    nm_total = alignment.NM + len_unaligned
    len_norm = alignment.blen + len_unaligned

    return nm_total / len_norm if len_norm > 0 else 1.0


@dataclass
class EdlibScoreResult:
    error: float
    cigar: str
    matched_target_sequence: str


def edlib_score(query_seq: str, target_seq: str, k: int = -1) -> Union[EdlibScoreResult, None]:
    """HW-align ``query_seq`` against a single ``target_seq`` and normalize to an
    error rate. Shared primitive for both assembly-window and read-based scoring;
    returns an ``EdlibScoreResult`` or ``None`` if no alignment was produced.

    ``k`` is edlib's max edit distance: alignments worse than ``k`` abort early
    and return ``None`` (edlib editDistance = -1). ``k=-1`` (default) is
    unbounded, preserving the assembly path's behavior; the read path passes a
    threshold-derived ``k`` so non-matching reads don't cost a full O(len*len)
    alignment -- the dominant cost when an SV is a read-mode miss."""
    result = edlib.align(query_seq.upper(), target_seq.upper(), mode="HW", task="path", k=k)
    if not result or result['editDistance'] < 0:
        return None

    edit_dist = result['editDistance']
    if result.get('locations') and result['locations'][0][0] is not None and result['locations'][0][1] is not None:
        loc = result['locations'][0]
        target_match_len = loc[1] - loc[0] + 1
        denominator = max(len(query_seq), target_match_len)
        matched_target_sequence = target_seq[loc[0]:loc[1] + 1]
    else:
        denominator = len(query_seq)
        matched_target_sequence = target_seq # TODO: check for correctness

    error_rate = edit_dist / denominator if denominator > 0 else 1.0
    return EdlibScoreResult(error=error_rate, cigar=result['cigar'], matched_target_sequence=matched_target_sequence)


def run_edlib_fallback(query_seq: str, chrom: str, location: int, initial_buffer: int, max_tolerance: int,
                       error_threshold: float, sample_dict: Dict[str, bytearray]) -> Union[EdlibScoreResult, None]:
    if sample_dict is None:
        return None

    target_keys = [k for k in sample_dict.keys() if k == chrom or k.startswith(f"{chrom}_")]
    if not target_keys:
        return None

    best_error = 1.0
    best_res = None
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

            res = edlib_score(query_seq, target_seq)
            if res is not None and res.error < best_error:
                best_error = res.error
                best_res = res

            if best_error <= error_threshold:
                return best_res

        if not expanded_any or current_tolerance >= max_tolerance:
            break

        current_tolerance = max(1, current_tolerance * 10)
        if current_tolerance > max_tolerance:
            current_tolerance = max_tolerance

    return best_res
