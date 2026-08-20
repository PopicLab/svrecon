"""Edlib-based validation: the expanding-window fallback search and the EdlibScorer."""
from typing import Dict, Union

from svrecon.constants import ValidationSource
from svrecon.reconstruct import Query
from svrecon.scorers.base import Scorer, QueryValidationInput
from svrecon.scorers.utils import Cigar, EdlibScoreResult, edlib_score, validate_segments_from_cigar
from svrecon.utils import load_fasta_to_bytes


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
