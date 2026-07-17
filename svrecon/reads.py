"""Read-based validation: thread-safe BAM reader and per-read edlib scoring."""
import logging
import os
import threading
from typing import Dict, List, Tuple, Union

import pysam

from svrecon.align import edlib_score
from svrecon.util import reverse_complement

logger = logging.getLogger(__name__)


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
