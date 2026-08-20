"""Edlib-based validation: HW alignment scoring, the expanding-window fallback search,
and the EdlibScorer. (`import edlib` below resolves to the external edlib package --
Python 3 imports are absolute by default.)"""
from dataclasses import dataclass
from typing import Dict, Union

import edlib

from svrecon.constants import ValidationSource
from svrecon.reconstruct import Query
from svrecon.scorers.base import Scorer, QueryValidationInput
from svrecon.scorers.cigar import Cigar, validate_segments_from_cigar
from svrecon.util import load_fasta_to_bytes


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


class EdlibScorer(Scorer):
    """Fallback for short queries mappy missed: expanding-window edlib search of a FASTA's bytes."""

    def __init__(self, config, chroms, fasta_path: str):
        super().__init__(config, chroms)
        self.min_edlib_query = config.min_edlib_query
        self.edlib_fallback_max_tolerance = config.edlib_fallback_max_tolerance
        self.sample_bytes = load_fasta_to_bytes(fasta_path, chroms)

    def score_query(self, query: Query) -> QueryValidationInput:
        if len(query.sequence) >= self.min_edlib_query:  # cost bound: the search is O(len * window)
            return QueryValidationInput(source=ValidationSource.EDLIB)

        lowest_pass_error = 1.0
        passed = False
        validating_seq = None
        segment_validation_results = []

        edlib_result = run_edlib_fallback(query.sequence, query.chrom, query.ref_start, query.buffer,
                                          self.edlib_fallback_max_tolerance, self.match_error_threshold,
                                          self.sample_bytes)
        if edlib_result is None:
            return QueryValidationInput(source=ValidationSource.EDLIB)

        if edlib_result.error <= self.match_error_threshold:
            segment_validation_results = validate_segments_from_cigar(
                Cigar.from_edlib(edlib_result.cigar), query.recon_segments,
                error_threshold=self.match_error_threshold)
            if all(s.passed for s in segment_validation_results):
                passed = True
                lowest_pass_error = edlib_result.error
                validating_seq = edlib_result.matched_target_sequence

        return QueryValidationInput(
            source=ValidationSource.EDLIB,
            passed=passed,
            lowest_pass_error=lowest_pass_error,
            lowest_error=edlib_result.error,
            validating_seq=validating_seq,
            segment_results=segment_validation_results,
        )
