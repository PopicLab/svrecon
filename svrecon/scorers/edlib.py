"""Edlib-based validation: the expanding-window fallback search and the EdlibScorer."""
from typing import Dict, List, Union

from svrecon.constants import SubseqReason, SubseqStatus, ValidationSource
from svrecon.reconstruct import Query
from svrecon.scorers.base import CigarValidationResult, Scorer, QueryValidation
from svrecon.scorers.utils import Cigar, EdlibScoreResult, edlib_score
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

    def score_query(self, query: Query) -> QueryValidation:
        if len(query) >= self.min_edlib_query:  # cost bound: the search is O(len * window)
            return QueryValidation(source=ValidationSource.EDLIB)

        edlib_result = run_edlib_fallback(query.sequence, query.chrom, query.ref_start, query.buffer,
                                          self.edlib_fallback_max_tolerance, self.match_error_threshold,
                                          self.sample_bytes)
        if edlib_result is None:  # nothing aligned at all
            return QueryValidation(source=ValidationSource.EDLIB)

        cigar = Cigar.from_edlib(edlib_result.cigar)
        cigar_results: List[CigarValidationResult] = self.score_cigar(cigar, query)
        passed = all(r.passed for r in cigar_results)

        if passed:  # a window aligned and every check passed
            status, reason = SubseqStatus.PASS, SubseqReason.PASS
        else:  # the best window aligned, but a check failed
            status, reason = SubseqStatus.MATCH, SubseqReason.CIGAR_FAILED

        return QueryValidation(
            source=ValidationSource.EDLIB,
            passed=passed,
            status=status,
            reason=reason,
            lowest_pass_error=edlib_result.error if passed else 1.0,
            lowest_error=edlib_result.error,
            best_matched_seq=edlib_result.matched_target_sequence,
            cigar_results=cigar_results,
            best_strand_match=1,  # edlib aligns forward only
            cigar=cigar,
        )
