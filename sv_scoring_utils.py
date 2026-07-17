import logging
import os
import threading
from pathlib import Path
from typing import Dict, List, Tuple, Union

import re
import edlib
import mappy
import pysam
import yaml
from pysam import VariantRecord

logger = logging.getLogger(__name__)

from typing import List, Tuple
import re


def edlib_to_cigartuples(cigar_str: str) -> List[Tuple[int, int]]:
    """
    Translates an Edlib character-based CIGAR string into BAM standard integer tuples.
    """
    op_map = {'=': 7, 'X': 8, 'I': 1, 'D': 2}
    return [(int(m.group(1)), op_map[m.group(2)]) for m in re.finditer(r'(\d+)([=XID])', cigar_str)]


def validate_junctions_from_cigar(cigartuples: List[Tuple[int, int]], junctions: List[int], q_st: int = 0,
                                  window: int = 150, error_threshold: float = 0.1) -> bool:
    """
    Calculates the local error rate within a specified window around structural variant junctions.
    Uses q_st to synchronize absolute junction indices with the relative CIGAR string.

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

        # Validate the error rate if the window contained sequence data
        if total_bases > 0:
            err_rate = errors / total_bases
            if err_rate > error_threshold:
                return False

    return True

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
        logger.debug(f'Full groovi config: {config_data}')
    return args


def get_start_stop(rec: VariantRecord) -> Tuple[int, int]:
    """Converts VCF start/stop coordinates to Python style"""
    start = rec.start - 1
    stop = rec.stop
    return start, stop


_RC_TRANS = str.maketrans('ACGTNacgtn', 'TGCANtgcan')


def reverse_complement(seq: Union[str, List[str]]) -> str:
    """Reverse complement of a DNA sequence (unknown bases -> N). Accepts a str or
    a list of single-character strings; always returns a str."""
    if not isinstance(seq, str):
        seq = ''.join(seq)
    return seq.translate(_RC_TRANS)[::-1]


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


class BamReader:
    """Process-wide, thread-safe BAM reader for read-based evaluation.

    A pysam.AlignmentFile handle is not safe for concurrent access, so a single
    shared handle is guarded by a lock: worker threads block on disk I/O exactly
    as needed (I/O is the bottleneck), while edlib alignment runs outside the
    lock in the caller. Requires a coordinate-sorted, indexed BAM (.bai); the
    index is built on construction if missing.
    """

    def __init__(self, bam_path: str):
        self.bam_path = bam_path
        if not (os.path.exists(bam_path + '.bai') or os.path.exists(bam_path + '.csi')):
            logger.info(f'No BAM index found; building one for {bam_path}')
            pysam.index(bam_path)
        self._bam = pysam.AlignmentFile(bam_path, 'rb')
        self._lock = threading.Lock()

    def candidate_read_seqs(self, chrom: str, start: int, end: int, flank: int,
                            max_reads: int) -> List[str]:
        """Full-molecule sequences of reads with ANY alignment overlapping
        ``[start - flank, end + flank]`` on ``chrom``.

        We deliberately cast a wide net rather than requiring a read to *span*
        the locus: a read that carries the SV is typically split into a primary
        + supplementary (chimeric) alignments, so no single alignment spans it.
        We therefore keep supplementary alignments, dedup by read name, and keep
        the longest sequence seen per read -- the primary alignment's
        soft-clipped record carries the complete molecule, while supplementary
        records are hard-clipped. The caller aligns the reconstruction against
        these full molecules with edlib and lets that judge support (we trust
        the aligner only to gather candidates, not to validate). Secondary and
        unmapped records are skipped (no usable / no full sequence).

        Collection stops once ``max_reads`` distinct reads have been gathered,
        which bounds work (and lock-hold time) on deep read pileups.

        Only the fetch + sequence copy is done under the lock; alignment is the
        caller's job."""
        lo = max(0, start - flank)
        hi = end + flank
        by_name = {}
        with self._lock:
            try:
                fetched = self._bam.fetch(chrom, lo, hi)
            except (ValueError, KeyError):
                return []  # chrom absent from BAM header
            for r in fetched:
                if r.is_unmapped or r.is_secondary:
                    continue
                seq = r.query_sequence
                if not seq:
                    continue
                name = r.query_name
                if name in by_name:
                    if len(seq) > len(by_name[name]):
                        by_name[name] = seq
                else:
                    by_name[name] = seq
                    if len(by_name) >= max_reads:
                        break
        return list(by_name.values())

    def close(self):
        with self._lock:
            self._bam.close()


def run_read_edlib(query_seq: str, read_seqs: List[str], error_threshold: float,
                   min_support: int = 1) -> Tuple[Union[Dict, None], int]:
    """Align ``query_seq`` against each candidate read with the shared
    ``edlib_score`` primitive.

    Returns ``(best_res, n_tried)``:
      * ``best_res`` -- the best ``{'error', 'cigar'}`` once ``min_support`` reads
        clear ``error_threshold`` (early exit), else the best within the 2x bound,
        else ``None``.
      * ``n_tried`` -- number of reads actually aligned (i.e. long enough to pass
        the length pre-filter). ``n_tried == 0`` means every candidate read was
        too short to host the alt -- distinct from "reads aligned but failed",
        so the caller can report it differently."""
    best_res = None
    best_error = 1.0
    support = 0
    n_tried = 0
    # A read can't contain the alt under HW alignment if it is shorter than the
    # alt by more than the error budget: the unmatched overhang alone forces
    # error >= (len(alt) - len(read)) / len(alt). Skip those reads.
    min_target_len = len(query_seq) * (1.0 - error_threshold)
    # Bound each edlib alignment so non-matching reads abort early instead of
    # computing a full O(len*len) matrix -- without this, a multi-kb alt that no
    # read supports grinds through every read at full cost. We allow up to 2x the
    # decision threshold so near-misses (threshold..2x) still compute and report
    # their error for diagnostics; only clearly-bad alignments (>2x) abort.
    k = max(1, int(2 * error_threshold * len(query_seq)))
    for target_seq in read_seqs:
        if len(target_seq) < min_target_len:
            continue
        n_tried += 1
        # A read molecule can be sequenced from either strand relative to the
        # reference-oriented alt, so try BOTH orientations and keep the better.
        candidates = [r for r in (edlib_score(query_seq, target_seq, k=k),
                                   edlib_score(query_seq, reverse_complement(target_seq), k=k))
                      if r is not None]
        if not candidates:
            continue
        res = min(candidates, key=lambda r: r['error'])
        if res['error'] < best_error:
            best_error = res['error']
            best_res = res
        if res['error'] <= error_threshold:
            support += 1
            if support >= min_support:
                return best_res, n_tried
    return best_res, n_tried
