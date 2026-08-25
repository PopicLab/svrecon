"""Assembly-based validation: per-chromosome mappy aligner cache and the AssemblyScorer."""
import hashlib
import logging
import os
from pathlib import Path
from typing import Dict, Iterable, List

import mappy
import pysam

from svrecon.constants import QueryValidationReason, QueryValidationStatus, ValidationSource
from svrecon.reconstruct import Query
from svrecon.scorers.base import Scorer, QueryValidation
from svrecon.scorers.utils import Cigar
from svrecon.utils import reverse_complement

logger = logging.getLogger(__name__)

# minimap2 --eqx: report =/X instead of M.
MM_F_EQX = 0x4000000

# Aligner build params, shared by the aligners
ALIGN_PARAMS = {
    'preset': 'map-hifi',
    'k': 15,
    'w': 5,
    'best_n': 100,
    'min_cnt': 1,
    'min_dp_score': 10,
    'min_chain_score': 1,
    'extra_flags': MM_F_EQX,
}


def get_chrom_aligner(sample_fasta: str, chrom: str, cache_dir: str, align_params: dict,
                      threads: int = 4) -> mappy.Aligner:
    """Loads a per-chromosome MMI index, building it first if it doesn't exist or the FASTA
    changed since it was built. Raises ValueError if the FASTA has no sequences for the chromosome."""
    mmi_path = os.path.join(cache_dir, f"{chrom}.mmi")
    # The cache is keyed on the FASTA's path, so a FASTA rewritten in place (a re-run
    # simulation, an updated assembly) would silently reuse an index of the old genome.
    source_path = os.path.join(cache_dir, f"{chrom}.source")
    fasta_stat = os.stat(sample_fasta)
    fasta_id = f'{fasta_stat.st_size},{fasta_stat.st_mtime_ns}'

    if os.path.exists(mmi_path) and os.path.exists(source_path):
        cached_id = open(source_path).read()
        if cached_id == fasta_id:
            logger.info(f"Loading existing index for {chrom} from {mmi_path} ({sample_fasta} {fasta_id})")
            return mappy.Aligner(mmi_path, **align_params)
        logger.warning(f"Stale index for {chrom}: built from {cached_id}, "
                       f"{sample_fasta} is now {fasta_id}. Rebuilding.")

    logger.info(f"Index not current for {chrom}. Extracting sequences...")
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
        raise ValueError(f"no sequences found for chromosome {chrom!r} in {sample_fasta}")

    logger.info(f"Building MMI index for {chrom} at {mmi_path}...")
    aligner = mappy.Aligner(chrom_fa_path, n_threads=threads, fn_idx_out=mmi_path, **align_params)
    with open(source_path, 'w') as f:
        f.write(fasta_id)

    if os.path.exists(chrom_fa_path):
        os.remove(chrom_fa_path)

    return aligner


def build_chrom_aligners(fasta_path: str, chroms: Iterable[str], cache_dir: Path) -> Dict[str, mappy.Aligner]:
    """Build/load per-chromosome mappy aligners for `fasta_path`, keyed by chromosome."""
    path_hash = hashlib.md5(str(Path(fasta_path).resolve()).encode('utf-8')).hexdigest()[:8]
    fa_cache = cache_dir / f'{Path(fasta_path).name}_{path_hash}'
    fa_cache.mkdir(parents=True, exist_ok=True)
    return {chrom: get_chrom_aligner(fasta_path, chrom, str(fa_cache), ALIGN_PARAMS, threads=32)
            for chrom in chroms}

class AssemblyScorer(Scorer):
    """Validates a query against a FASTA (sample assembly, or the reference for the
    ambiguity check) via per-chromosome mappy alignments."""

    def __init__(self, config, chroms, fasta_path: str, forward_match_only: bool = False):
        super().__init__(config, chroms)
        self.forward_match_only = forward_match_only
        self.aligners = build_chrom_aligners(fasta_path, chroms, config.cache_dir)
        logger.info(f'Initialized per-chromosome aligners for {fasta_path}')

    def score_query(self, query: Query) -> QueryValidation:
        lowest_pass_error = 1.0
        lowest_error = 1.0
        passed = False
        inconclusive = False
        best_cigar_results = []
        best_strand_match = None
        best_cigar = None
        best_matched_seq = None

        chrom_aligner = self.aligners[query.chrom]

        # gather and filter candidate alignments
        alignments: List[mappy.Alignment] = list(chrom_aligner.map(query.sequence)) # .map returns possible alignments, esp repetitive alignments?
        alignments = [a for a in alignments if abs(a.r_st - query.ref_start) <= self.location_tolerance]
        if self.forward_match_only:
            alignments = [a for a in alignments if a.strand == 1]

        # pick best candidate alignment by bulk error; every configured check must pass # TODO: ensure mapq 60. something thats repetitive ins unknown
        for alignment in alignments:
            cigar = Cigar.from_mappy(alignment, len(query))
            cigar_status, cigar_results = self.score_cigar(cigar, query)
            match_err = cigar.error_rate
            lowest_error = min(lowest_error, match_err)

            matched_seq = chrom_aligner.seq(alignment.ctg, alignment.r_st, alignment.r_en)
            if alignment.strand == -1:
                matched_seq = reverse_complement(matched_seq)  # orient to the forward query
            if cigar_status is QueryValidationStatus.INCONCLUSIVE:
                inconclusive = True
            if cigar_status is QueryValidationStatus.PASS:
                passed = True
                if match_err < lowest_pass_error:
                    lowest_pass_error = match_err
                    best_strand_match = alignment.strand
                    best_cigar_results = cigar_results
                    best_cigar = cigar
                    best_matched_seq = matched_seq
            if not passed and match_err <= lowest_error:
                best_cigar_results = cigar_results
                best_cigar = cigar
                best_matched_seq = matched_seq

        if passed:
            status, reason = QueryValidationStatus.PASS, QueryValidationReason.PASS
        elif inconclusive:  # aligned, but a check could not judge the query
            status, reason = QueryValidationStatus.INCONCLUSIVE, QueryValidationReason.OTHER
        elif best_cigar_results:  # aligned to the sample, but a check failed
            status, reason = QueryValidationStatus.FAIL, QueryValidationReason.CIGAR_FAILED
        else:  # nothing aligned at all
            status, reason = QueryValidationStatus.FAIL, QueryValidationReason.OTHER

        return QueryValidation(
            source=ValidationSource.ASSEMBLY,
            passed=passed,
            status=status,
            reason=reason,
            lowest_pass_error=lowest_pass_error,
            lowest_error=lowest_error,
            best_matched_seq=best_matched_seq,
            cigar_results=best_cigar_results,
            best_strand_match=best_strand_match,
            cigar=best_cigar,
        )
