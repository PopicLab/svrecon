"""Read-based validation: thread-safe BAM reader, per-read edlib scoring, and the ReadScorer."""
import logging
import os
import threading
from typing import List, Optional
import pysam

from svrecon.constants import QueryValidationReason, QueryValidationStatus, ValidationSource
from svrecon.reconstruct import Query
from svrecon.scorers.base import Scorer, QueryValidation
from svrecon.scorers.utils import Cigar, EdlibScoreResult, edlib_score
from svrecon.utils import reverse_complement

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

    def candidate_read_seqs(self, chrom: str, start: int, end: int, max_reads: int) -> List[str]:
        """Fetch full-molecule sequences of reads with ANY alignment overlapping
        ``[start, end]`` on ``chrom``."""
        lo = max(0, start)
        hi = end
        by_name = {}
        with self._lock:
            fetched = self._bam.fetch(chrom, lo, hi)
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


def run_read_edlib(query_seq: str, read_seqs: List[str], error_threshold: float) -> List[EdlibScoreResult]:
    """Return the EdlibScoreResult for each read (forward and RC tried) that contains
    ``query_seq`` with at most ``error_threshold`` error."""
    passing: List[EdlibScoreResult] = []
    # A read can't contain the alt under HW alignment if it is shorter than the
    # alt by more than the error budget: the unmatched overhang alone forces
    # error >= (len(alt) - len(read)) / len(alt). Skip those reads.
    min_target_len = len(query_seq) * (1.0 - error_threshold)
    # Bound each edlib alignment so non-matching reads abort early instead of
    # computing a full O(len*len) matrix -- without this, a multi-kb alt that no
    # read supports grinds through every read at full cost.
    k = max(1, int(2 * error_threshold * len(query_seq)))
    for target_seq in read_seqs:
        if len(target_seq) < min_target_len:
            continue
        # A read molecule can be sequenced from either strand relative to the
        # reference-oriented alt, so try BOTH orientations and keep the better.
        candidates = [r for r in (edlib_score(query_seq, target_seq, k=k),
                                  edlib_score(query_seq, reverse_complement(target_seq), k=k))
                      if r is not None]
        if not candidates:
            continue
        res = min(candidates, key=lambda r: r.error)
        if res.error <= error_threshold:
            passing.append(res)
    return passing


class ReadScorer(Scorer):
    """Validates a query against real reads: does any single read contain the allele?"""

    def __init__(self, config, chroms):
        super().__init__(config, chroms, error_threshold=config.read_error_threshold)
        self.read_error_threshold = config.read_error_threshold
        self.max_reads_per_site = config.max_reads_per_site
        self.bam_reader = BamReader(config.bam)

    def score_query(self, query: Query) -> QueryValidation:
        # collect candidate reads that overlap ref region of interest, then filter by size
        breakpoints = [pos for start, end in query.ref_segments for pos in (start, end)]
        bp_start, bp_stop = min(breakpoints), max(breakpoints)
        reads = self.bam_reader.candidate_read_seqs(
            query.chrom, bp_start, bp_stop, max_reads=self.max_reads_per_site)
        reads = [r for r in reads if len(r) >= len(query)]
        if not reads:  # no read long enough to span the allele -> untestable, not contradicted
            return QueryValidation(source=ValidationSource.READS,
                                        status=QueryValidationStatus.INCONCLUSIVE,
                                        reason=QueryValidationReason.OTHER)

        passed = False
        inconclusive = False
        lowest_pass_error = 1.0
        lowest_error = 1.0
        best_cigar_results = []  # the adopted read's checks on a pass, else the best-scoring read's
        best_cigar: Optional[Cigar] = None
        best_matched_seq: Optional[str] = None

        read_results = run_read_edlib(query.sequence, reads, self.read_error_threshold)
        for read_res in read_results:
            lowest_error = min(lowest_error, read_res.error)
            cigar: Cigar = read_res.cigar
            cigar_status, cigar_results = self.score_cigar(cigar, query)
            if cigar_status is QueryValidationStatus.INCONCLUSIVE:
                inconclusive = True
            if cigar_status is QueryValidationStatus.PASS and read_res.error < lowest_pass_error:
                passed = True
                lowest_pass_error = read_res.error
                best_cigar_results = cigar_results
                best_cigar = cigar
                best_matched_seq = read_res.matched_target_sequence
            elif not passed and read_res.error <= lowest_error:
                best_cigar_results = cigar_results
                best_cigar = cigar
                best_matched_seq = read_res.matched_target_sequence  # kept so an aligned fail can still be inspected

        if passed:
            status, reason = QueryValidationStatus.PASS, QueryValidationReason.PASS
        elif inconclusive:  # a read aligned, but a check could not judge the query
            status, reason = QueryValidationStatus.INCONCLUSIVE, QueryValidationReason.OTHER
        elif read_results:  # a read aligned, but none passed every check
            status, reason = QueryValidationStatus.FAIL, QueryValidationReason.CIGAR_FAILED
        else:  # no read within the error budget
            status, reason = QueryValidationStatus.FAIL, QueryValidationReason.OTHER

        return QueryValidation(source=ValidationSource.READS, passed=passed,
                                    status=status, reason=reason,
                                    lowest_pass_error=lowest_pass_error, lowest_error=lowest_error,
                                    best_matched_seq=best_matched_seq, cigar_results=best_cigar_results,
                                    cigar=best_cigar)
