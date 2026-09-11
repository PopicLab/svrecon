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

    def candidate_read_seqs(self, chrom: str, start: int, end: int, max_reads: int) -> List[pysam.AlignedSegment]:
        """Fetch full-molecule sequences of reads with ANY alignment overlapping
        ``[start, end]`` on ``chrom``."""
        lo = max(0, start)
        hi = end
        by_name = {}
        with self._lock:
            fetched = self._bam.fetch(chrom, lo, hi)
            for read in fetched:
                if read.is_unmapped or read.is_secondary:
                    continue
                seq = read.query_sequence
                if not seq:
                    continue
                name = read.query_name
                if name in by_name:
                    if len(seq) > len(by_name[name]):
                        by_name[name] = read
                else:
                    by_name[name] = read
                    if len(by_name) >= max_reads:
                        break
        return list(by_name.values())


def run_read_edlib(query_seq: str, reads: List[pysam.AlignedSegment], error_threshold: float) -> List[EdlibScoreResult]:
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
    read_seqs = [read.query_sequence for read in reads]
    for target_seq in read_seqs:
        if len(target_seq) < min_target_len:
            continue
        # A read molecule can be sequenced from either strand relative to the
        # reference-oriented alt, so try BOTH orientations and keep the better.
        candidates = [edlib_result for edlib_result in (edlib_score(query_seq, target_seq, k=k),
                                  edlib_score(query_seq, reverse_complement(target_seq), k=k))
                      if edlib_result is not None]
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
        reads = [r for r in reads if r.query_length >= len(query)]

        # no read long enough to span the allele -> inconclusive
        if not reads:  
            return QueryValidation(source=ValidationSource.READS,
                                        status=QueryValidationStatus.INCONCLUSIVE,
                                        reason=QueryValidationReason.NO_SPANNING_READS)
        # filter by read error threshold
        filtered_read_results = run_read_edlib(query.sequence, reads, self.read_error_threshold)

        # if no reads passes the error threshold, score and fail against first read
        if not filtered_read_results:
            read = reads[0]
            query_positions = [
                q_pos for q_pos, ref_pos in read.get_aligned_pairs(matches_only=True)
                if bp_start <= ref_pos < bp_stop
            ]
            if query_positions:
                block_sequence = read.query_sequence[query_positions[0]:query_positions[-1] + 1]
                edlib_results = edlib_score(query.sequence, block_sequence, mode='NW')
                cigar_status, cigar_results = self.score_cigar(edlib_results.cigar, query)

                return QueryValidation(source=ValidationSource.READS, passed=False,
                                    status=cigar_status, reason=QueryValidationReason.CIGAR_FAILED if cigar_status == QueryValidationStatus.FAIL else QueryValidationReason.CIGAR_INCONCLUSIVE,
                                    lowest_pass_error=1, lowest_error=edlib_results.cigar.whole_query_error_rate,
                                    best_matched_seq=block_sequence, cigar_results=cigar_results,
                                    cigar=edlib_results.cigar)
            else:
                # no matched base in [bp_start, bp_stop) 
                return QueryValidation(source=ValidationSource.READS,
                                       status=QueryValidationStatus.FAIL,
                                       reason=QueryValidationReason.CIGAR_FAILED)

        # score candidate read cigars
        cigar_scoring_results = [self.score_cigar(read_result.cigar, query) for read_result in filtered_read_results]

        # if all cigars are inconclusive, return inconclusive
        if all(cigar_status is QueryValidationStatus.INCONCLUSIVE for cigar_status, _ in cigar_scoring_results):
            return QueryValidation(source=ValidationSource.READS,
                                   status=QueryValidationStatus.INCONCLUSIVE,
                                   reason=QueryValidationReason.CIGAR_INCONCLUSIVE)
        
        # score remaining candidate reads
        passed = False
        lowest_pass_error = 1.0
        lowest_error = 1.0
        best_cigar_results = []
        best_cigar: Optional[Cigar] = None
        best_matched_seq: Optional[str] = None

        for i, read_res in enumerate(filtered_read_results):  
            cigar_status, cigar_results = cigar_scoring_results[i]
            if cigar_status == QueryValidationStatus.PASS:
                passed = True
                if read_res.error <= lowest_pass_error:
                    lowest_pass_error = read_res.error
                    best_cigar = read_res.cigar
                    best_cigar_results = cigar_results
                    best_matched_seq = read_res.matched_target_sequence
                    lowest_error = min(lowest_error, lowest_pass_error)
            elif not passed and read_res.error <= lowest_error:
                lowest_error = read_res.error
                best_cigar = read_res.cigar
                best_cigar_results = cigar_results
                best_matched_seq = read_res.matched_target_sequence

        if passed:
            status, reason = QueryValidationStatus.PASS, QueryValidationReason.PASS
        else: 
            status, reason = QueryValidationStatus.FAIL, QueryValidationReason.CIGAR_FAILED

        return QueryValidation(source=ValidationSource.READS, passed=passed,
                                    status=status, reason=reason,
                                    lowest_pass_error=lowest_pass_error, lowest_error=lowest_error,
                                    best_matched_seq=best_matched_seq, cigar_results=best_cigar_results,
                                    cigar=best_cigar)
