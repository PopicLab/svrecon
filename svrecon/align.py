"""Alignment + error-rate primitives: mappy aligner caching, CIGAR/junction validation, edlib scoring."""
import logging
import os
import re
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


def validate_junctions_from_cigar(cigartuples: List[Tuple[int, int]], junctions: List[int], q_st: int = 0,
                                  window: int = 150, error_threshold: float = 0.1) -> List[Dict]:
    """
    Calculates the local error rate within a specified window around structural variant junctions.
    Uses q_st to synchronize absolute junction indices with the relative CIGAR string.

    Returns an ordered list of per-junction results, one ``{'error', 'passed'}`` dict per
    in-scope junction in the order checked. Evaluation short-circuits on the first junction
    whose local error exceeds ``error_threshold``: that failing junction is the last entry and
    later in-scope junctions are not checked. A junction whose window contains no assessable
    bases is recorded as ``{'error': None, 'passed': True}``. Callers derive the overall verdict
    as ``all(j['passed'] for j in results)`` (an empty list -- no in-scope junctions -- passes).

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

    # Conceptually group the operations
    # Operations that consume space in the simulated query sequence
    QUERY_CONSUMING_OPS = {MATCH, INS, SOFT_CLIP, SEQ_MATCH, SEQ_MISMATCH}

    # Operations that only represent gaps in the target assembly
    TARGET_ONLY_OPS = {DEL, REF_SKIP}

    # Calculate total query length consumed by this specific CIGAR
    cigar_q_len = sum(length for length, op in cigartuples if op in QUERY_CONSUMING_OPS)
    q_en = q_st + cigar_q_len

    results = []
    for j_idx in junctions:
        # Only validate junctions that fall within the scope of this alignment segment.
        # This prevents "False Misses" when Mappy splits chimeric alignments.
        if j_idx < q_st - window or j_idx > q_en + window:
            continue

        # Translate the absolute query junction index to a relative position within this CIGAR
        relative_j_idx = j_idx - q_st

        window_start = max(0, relative_j_idx - window)
        window_end = min(cigar_q_len, relative_j_idx + window)

        query_cursor = 0
        errors = 0
        total_bases = 0

        for length, op in cigartuples:
            # Stop processing if the cursor has moved past the evaluation window
            if query_cursor > window_end and op not in TARGET_ONLY_OPS:
                break

            # --- TARGET-ONLY OPERATIONS (Deletions / Skips) ---
            if op in TARGET_ONLY_OPS:
                # If a deletion occurs while the cursor is inside the window, it is recorded as an error.
                if window_start <= query_cursor < window_end:
                    errors += length
                    total_bases += length
                continue

            # --- QUERY-CONSUMING OPERATIONS ---
            op_start = query_cursor
            op_end = query_cursor + length

            # Calculate the overlap between this CIGAR operation and the evaluation window
            overlap_start = max(window_start, op_start)
            overlap_end = min(window_end, op_end)
            overlap_len = max(0, overlap_end - overlap_start)

            if overlap_len > 0:
                if op in (MATCH, SEQ_MATCH):
                    total_bases += overlap_len
                elif op in (INS, SOFT_CLIP, SEQ_MISMATCH):
                    errors += overlap_len
                    total_bases += overlap_len

            # Advance the query cursor
            if op in QUERY_CONSUMING_OPS:
                query_cursor += length

        # Record the per-junction result. Windows with no assessable bases cannot be
        # judged, so they are treated as passing (matching the original skip behavior).
        if total_bases > 0:
            err_rate = errors / total_bases
            passed = err_rate <= error_threshold
            # Store as a native JSON number at 4 sig figs; small rates keep precision and
            # serialize in exponent form (e.g. 8.6e-06) rather than flattening to 0.
            results.append({'error': float(f'{err_rate:.4g}'), 'passed': passed})
            # Short-circuit on the first failing junction (original return-False behavior).
            if not passed:
                break
        else:
            results.append({'error': None, 'passed': True})

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
