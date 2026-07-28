"""Alignment + error-rate primitives: mappy aligner caching, CIGAR/junction validation, edlib scoring."""
import logging
import os
import re
from dataclasses import dataclass
from typing import Dict, List, Tuple, Union

import edlib
import mappy
import pysam

logger = logging.getLogger(__name__)


def edlib_to_cigartuples(cigar_str: str) -> List[Tuple[int, int]]:
    """
    Translates an Edlib character-based CIGAR string into BAM standard integer tuples.
    """
    op_map = {'=': 7, 'X': 8, 'I': 1, 'D': 2}
    return [(int(m.group(1)), op_map[m.group(2)]) for m in re.finditer(r'(\d+)([=XID])', cigar_str)]


@dataclass
class SeqJunctionsValidationResult:
    error: float
    passed: bool

    def jsonify(self) -> dict:
        return {'error': self.error, 'passed': self.passed}


def validate_junctions_from_cigar(cigartuples: List[Tuple[int, int]], junctions: List[int], q_st: int = 0,
                                  q_en: int = None, query_len: int = None, strand: int = 1,
                                  window: int = 150, error_threshold: float = 0.1) -> List[SeqJunctionsValidationResult]:
    """
    Calculates the local error rate within a window around each structural-variant junction.

    The window is defined in FORWARD-QUERY coordinates ``[j - window, j + window]`` (clamped to
    ``[0, query_len]``). Two subtleties the aligner introduces are handled explicitly:

    * Soft-clipped flanks. mappy's CIGAR covers only the aligned region ``[q_st, q_en)``; any
      query base in the window that the aligner clipped (outside that region) is a base where the
      junction adjacency is NOT actually spanned, so it is counted as an error. Without this a
      novel junction sitting at the clip boundary (e.g. dupINVdup's outer ``A|c`` / ``a|C``, where
      only the inverted core matches the reference) would read as clean 0.0 error.
    * Reverse-strand orientation. A reverse hit's CIGAR walks ``RC(query)``, so forward-query
      offset ``p`` sits at CIGAR-query cursor ``q_en - 1 - p``. The forward slice of the window is
      therefore mapped to the mirrored cursor slice before the CIGAR is crawled.

    ``q_en`` / ``query_len`` / ``strand`` default to the forward, fully-consumed case
    (``q_en = q_st + cigar_q_len``, ``query_len = q_en``, ``strand = 1``) so edlib callers -- which
    consume the whole query on the forward strand -- can omit them.

    Returns an ordered list of per-junction results, one ``SeqJunctionsValidationResult`` per
    in-scope junction in the order checked. Evaluation short-circuits on the first junction
    whose local error exceeds ``error_threshold``: that failing junction is the last entry and
    later in-scope junctions are not checked. A junction whose window collapses to no bases is
    recorded as ``SeqJunctionsValidationResult(error=0.0, passed=True)``. Callers derive the
    overall verdict as ``all(j.passed for j in results)`` (an empty list -- no in-scope
    junctions -- passes).

    Official SAM/BAM CIGAR Specification: https://samtools.github.io/hts-specs/SAMv1.pdf
    """
    # Define CIGAR operation constants for readability
    MATCH = 0  # M: Alignment match (can be sequence match or mismatch)
    INS = 1  # I: Insertion to the reference
    DEL = 2  # D: Deletion from the reference
    REF_SKIP = 3  # N: Skipped region from the reference
    SOFT_CLIP = 4  # S: Soft clipping (clipped sequences present in query)
    HARD_CLIP = 5  # H: Hard clipping (clipped sequences NOT present in query)
    PAD = 6  # P: Padding (silent deletion from padded reference)
    SEQ_MATCH = 7  # =: Exact sequence match
    SEQ_MISMATCH = 8  # X: Exact sequence mismatch

    # Operations that consume space in the simulated query sequence
    QUERY_CONSUMING_OPS = {MATCH, INS, SOFT_CLIP, SEQ_MATCH, SEQ_MISMATCH}
    # Operations that only represent gaps in the target assembly
    TARGET_ONLY_OPS = {DEL, REF_SKIP}
    # Query-consuming operations that constitute local error
    ERROR_OPS = {INS, SOFT_CLIP, SEQ_MISMATCH}

    # Calculate total query length consumed by this specific CIGAR
    cigar_q_len = sum(length for length, op in cigartuples if op in QUERY_CONSUMING_OPS)
    if q_en is None:
        q_en = q_st + cigar_q_len
    if query_len is None:
        query_len = q_en

    results = []
    for j_idx in junctions:
        # Only validate junctions that fall within the scope of this alignment segment.
        # This prevents "False Misses" when Mappy splits chimeric alignments.
        if j_idx < q_st - window or j_idx > q_en + window:
            continue

        # Window in FORWARD-QUERY coordinates, clamped to the real query.
        w_lo = max(0, j_idx - window)
        w_hi = min(query_len, j_idx + window)
        if w_hi <= w_lo:
            results.append(SeqJunctionsValidationResult(error=0.0, passed=True))
            continue

        total_bases = w_hi - w_lo
        # Clipped query bases inside the window: the adjacency is not spanned there -> error.
        errors = max(0, min(w_hi, q_st) - w_lo) + max(0, w_hi - max(w_lo, q_en))

        # Aligned slice of the window, mapped from forward-query into CIGAR-cursor coordinates.
        a_lo = max(w_lo, q_st)
        a_hi = min(w_hi, q_en)
        if a_hi > a_lo:
            if strand == -1:
                # Reverse hit: CIGAR walks RC(query); mirror the forward slice.
                c_lo = q_en - a_hi
                c_hi = q_en - a_lo
            else:
                c_lo = a_lo - q_st
                c_hi = a_hi - q_st

            query_cursor = 0
            for length, op in cigartuples:
                # Stop once the cursor has moved past the evaluation window.
                if query_cursor >= c_hi and op not in TARGET_ONLY_OPS:
                    break

                if op in TARGET_ONLY_OPS:
                    # A deletion inside the window is recorded as error.
                    if c_lo <= query_cursor < c_hi:
                        errors += length
                        total_bases += length
                    continue

                if op in QUERY_CONSUMING_OPS:
                    if op in ERROR_OPS:
                        overlap = max(0, min(c_hi, query_cursor + length) - max(c_lo, query_cursor))
                        errors += overlap
                    query_cursor += length

        err_rate = errors / total_bases
        passed = err_rate <= error_threshold
        # Store as a native JSON number at 4 sig figs; small rates keep precision and
        # serialize in exponent form (e.g. 8.6e-06) rather than flattening to 0.
        results.append(SeqJunctionsValidationResult(error=float(f'{err_rate:.4g}'), passed=passed))
        # Short-circuit on the first failing junction (original return-False behavior).
        if not passed:
            break

    return results


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


def check_match(alignment, query):
    """Calculate the independent error rate of a mappy alignment with the query sequence"""
    aligned_query_segment_length = alignment.q_en - alignment.q_st
    len_unaligned = len(query) - aligned_query_segment_length

    nm_total = alignment.NM + len_unaligned
    len_norm = alignment.blen + len_unaligned

    return nm_total / len_norm if len_norm > 0 else 1.0


def edlib_score(query_seq: str, target_seq: str, k: int = -1) -> Union[Dict, None]:
    """HW-align ``query_seq`` against a single ``target_seq`` and normalize to an
    error rate. Shared primitive for both assembly-window and read-based scoring;
    returns ``{'error', 'cigar'}`` or ``None`` if no alignment was produced.

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
    else:
        denominator = len(query_seq)

    error_rate = edit_dist / denominator if denominator > 0 else 1.0
    return {'error': error_rate, 'cigar': result['cigar']}


def run_edlib_fallback(query_seq: str, chrom: str, location: int, initial_buffer: int, max_tolerance: int,
                       error_threshold: float, sample_dict: Dict[str, bytearray]) -> Union[Dict, None]:
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
            if res is not None and res['error'] < best_error:
                best_error = res['error']
                best_res = res

            if best_error <= error_threshold:
                return best_res

        if not expanded_any or current_tolerance >= max_tolerance:
            break

        current_tolerance = max(1, current_tolerance * 10)
        if current_tolerance > max_tolerance:
            current_tolerance = max_tolerance

    return best_res
