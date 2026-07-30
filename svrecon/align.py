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

# CIGAR operation codes (SAM/BAM spec), in code order.
_MATCH, _INS, _DEL, _REF_SKIP, _SOFT_CLIP, _HARD_CLIP, _PAD, _SEQ_MATCH, _SEQ_MISMATCH = range(9)
# Op groupings used when scanning a CIGAR window.
_QUERY_CONSUMING_OPS = {_MATCH, _INS, _SOFT_CLIP, _SEQ_MATCH, _SEQ_MISMATCH}  # advance the query cursor
_TARGET_ONLY_OPS = {_DEL, _REF_SKIP}                                          # query gaps (deletions)
_ERROR_OPS = {_INS, _SOFT_CLIP, _SEQ_MISMATCH}                                # query bases that aren't clean matches


def edlib_to_cigartuples(cigar_str: str) -> List[Tuple[int, int]]:
    """
    Translates an Edlib character-based CIGAR string into BAM standard integer tuples.
    """
    op_map = {'=': 7, 'X': 8, 'I': 1, 'D': 2}
    return [(int(m.group(1)), op_map[m.group(2)]) for m in re.finditer(r'(\d+)([=XID])', cigar_str)]


@dataclass
class SeqSegmentValidationResult:
    error: float
    passed: bool

    def jsonify(self) -> dict:
        return {'error': self.error, 'passed': self.passed}


def validate_segments_from_cigar(cigartuples: List[Tuple[int, int]], segments: List[Tuple[int, int]], q_st: int = 0,
                                  q_en: int = None, query_len: int = None, strand: int = 1,
                                  radius: int = 150, error_threshold: float = 0.1) -> List[SeqSegmentValidationResult]:
    """
    Using the cigar summary of an alignment, verify the error rate for each of the provided segments is within tolerance.

    Each segment is a forward-query range ``[start, end)`` checked over the window
    ``[start - radius, end + radius]`` (clamped to ``[0, query_len]``). Any length works, so a pure
    deletion, which has no query extent, is just a zero-length segment ``(i, i)`` whose window is
    the ``radius`` context around its join. Two aligner subtleties are handled explicitly:

    * Soft-clipped flanks: query bases in the window outside the aligned region ``[q_st, q_en)`` are
      not spanned by the alignment and count as error (so a segment sitting at a clip boundary does
      not read as clean 0.0).
    * Reverse-strand orientation: a reverse hit's CIGAR walks ``RC(query)``, so the forward-query
      window is mirrored to CIGAR-cursor coordinates before the CIGAR is crawled.

    ``q_en`` / ``query_len`` / ``strand`` default to the forward, fully-consumed case so edlib
    callers (which consume the whole query on the forward strand) can omit them.

    :param cigartuples: BAM-style ``(length, op)`` tuples for the alignment.
    :param segments: List of ``(start, end)`` forward-query ranges to validate. A deletion is a
        zero-length ``(i, i)`` range.
    :param radius: context buffer window extended left and right of the segment range.

    Returns an ordered list of :class:`SeqSegmentValidationResult`, one per in-scope segment,
    short-circuiting on the first segment whose local error exceeds ``error_threshold``. A segment
    whose window collapses to no bases is recorded as ``error=0.0, passed=True``. Callers derive the
    overall verdict as ``all(s.passed for s in results)`` (an empty list passes).

    Official SAM/BAM CIGAR Specification: https://samtools.github.io/hts-specs/SAMv1.pdf
    """
    cigar_query_length = sum(length for length, op in cigartuples if op in _QUERY_CONSUMING_OPS)
    if q_en is None:
        q_en = q_st + cigar_query_length
    if query_len is None:
        query_len = q_en

    results = []
    for seg_start, seg_end in segments:
        # Skip segments outside this alignment's scope (a chimeric split can land one elsewhere).
        if seg_end < q_st - radius or seg_start > q_en + radius:
            continue

        # The segment plus `radius` of context, in forward-query coordinates, clamped to the query.
        window_start = max(0, seg_start - radius)
        window_end = min(query_len, seg_end + radius)
        if window_end <= window_start:
            results.append(SeqSegmentValidationResult(error=0.0, passed=True))
            continue

        error_rate = _window_error_rate(cigartuples, window_start, window_end, q_st, q_en, strand)
        passed = error_rate <= error_threshold
        results.append(SeqSegmentValidationResult(error=float(f'{error_rate:.4g}'), passed=passed))
        if not passed:  # short-circuit on the first failing segment
            break

    return results


def _window_error_rate(cigartuples: List[Tuple[int, int]], window_start: int, window_end: int,
                       q_st: int, q_en: int, strand: int) -> float:
    """Fraction of the forward-query window ``[window_start, window_end)`` that is not a clean
    aligned match. Two sources of error: query bases the aligner clipped (outside ``[q_st, q_en)``),
    and mismatch / insertion / deletion bases the CIGAR reports inside the window. A reverse hit's
    CIGAR walks ``RC(query)``, so the window is mirrored into CIGAR-cursor coordinates before it is
    scanned."""
    window_bases = window_end - window_start
    # Bases in the window the alignment never spanned (soft-clipped flanks).
    error_bases = max(0, min(window_end, q_st) - window_start) + max(0, window_end - max(window_start, q_en))

    aligned_start = max(window_start, q_st)
    aligned_end = min(window_end, q_en)
    if aligned_end > aligned_start:
        if strand == -1:
            cursor_start, cursor_end = q_en - aligned_end, q_en - aligned_start
        else:
            cursor_start, cursor_end = aligned_start - q_st, aligned_end - q_st

        query_cursor = 0
        for length, op in cigartuples:
            if query_cursor >= cursor_end and op not in _TARGET_ONLY_OPS:
                break  # past the window

            if op in _TARGET_ONLY_OPS:
                # A deletion inside the window is error and widens the assessed span.
                if cursor_start <= query_cursor < cursor_end:
                    error_bases += length
                    window_bases += length
                continue

            if op in _QUERY_CONSUMING_OPS:
                if op in _ERROR_OPS:
                    error_bases += max(0, min(cursor_end, query_cursor + length) - max(cursor_start, query_cursor))
                query_cursor += length

    return error_bases / window_bases


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
