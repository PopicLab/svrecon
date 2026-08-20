"""Assembly-based validation: per-chromosome mappy aligner cache and the AssemblyScorer."""
import hashlib
import logging
import os
from pathlib import Path
from typing import Dict, Iterable, List

import mappy
import pysam

from svrecon.constants import ValidationSource
from svrecon.reconstruct import Query
from svrecon.scorers.base import Scorer, QueryValidationInput
from svrecon.scorers.cigar import Cigar, SegmentValidation, validate_segments_from_cigar
from svrecon.util import reverse_complement

logger = logging.getLogger(__name__)

# Aligner build params, shared by the assembly aligners and (when
# --check_reference is on) the reference aligners.
ALIGN_PARAMS = {
    'preset': 'map-hifi',
    'k': 15,
    'w': 5,
    'best_n': 100,
    'min_cnt': 1,
    'min_dp_score': 10,
    'min_chain_score': 1,
}


def get_chrom_aligner(sample_fasta: str, chrom: str, cache_dir: str, align_params: dict,
                      threads: int = 4) -> mappy.Aligner:
    """Loads a per-chromosome MMI index, building it first if it doesn't exist.
    Raises ValueError if the FASTA has no sequences for the chromosome."""
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
        raise ValueError(f"no sequences found for chromosome {chrom!r} in {sample_fasta}")

    logger.info(f"Building MMI index for {chrom} at {mmi_path}...")
    aligner = mappy.Aligner(chrom_fa_path, n_threads=threads, fn_idx_out=mmi_path, **align_params)

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


# TODO: rename to something like get normalized score, change query param from str to query object
def check_match(alignment, query):
    """Calculate the independent error rate of a mappy alignment with the query sequence"""
    aligned_query_segment_length = alignment.q_en - alignment.q_st
    len_unaligned = len(query) - aligned_query_segment_length

    nm_total = alignment.NM + len_unaligned
    len_norm = alignment.blen + len_unaligned

    return nm_total / len_norm if len_norm > 0 else 1.0


class AssemblyScorer(Scorer):
    """Validates a query against a FASTA (sample assembly, or the reference for the
    ambiguity check) via per-chromosome mappy alignments."""

    def __init__(self, config, chroms, fasta_path: str, forward_match_only: bool = False):
        super().__init__(config, chroms)
        self.forward_match_only = forward_match_only
        self.aligners = build_chrom_aligners(fasta_path, chroms, config.cache_dir)
        logger.info(f'Initialized per-chromosome aligners for {fasta_path}')

    def score_query(self, query: Query) -> QueryValidationInput:
        lowest_pass_error = 1.0
        lowest_error = 1.0
        passed = False
        validating_seq = None
        best_segment_validation_results = []
        best_strand_match = None

        aligner = self.aligners.get(query.chrom)
        if aligner is None:
            return QueryValidationInput(source=ValidationSource.ASSEMBLY)

        # gather and filter candidate alignments
        alignments: List[mappy.Alignment] = list(aligner.map(query.sequence)) # .map returns possible alignments, esp repetitive alignments?
        alignments = [a for a in alignments if abs(a.r_st - query.ref_start) <= self.location_tolerance]
        if self.forward_match_only:
            alignments = [a for a in alignments if a.strand == 1]

        # pick best candidate alignment based on overall and segment normalized match error
        for alignment in alignments:
            match_err = check_match(alignment, query.sequence)
            lowest_error = min(lowest_error, match_err)
            if match_err > self.match_error_threshold:
                continue

            segment_validation_results: List[SegmentValidation] = \
                validate_segments_from_cigar(Cigar.from_mappy(alignment, len(query.sequence)),
                                             query.recon_segments,
                                             error_threshold=self.match_error_threshold)
            if all(s.passed for s in segment_validation_results):
                passed = True
                if match_err < lowest_pass_error:
                    lowest_pass_error = match_err
                    best_strand_match = alignment.strand   # record the best passing match's strand
                    best_segment_validation_results = segment_validation_results
                    matched = aligner.seq(alignment.ctg, alignment.r_st, alignment.r_en)
                    validating_seq = reverse_complement(matched) if alignment.strand == -1 else matched
            if not passed and match_err <= lowest_error:
                best_segment_validation_results = segment_validation_results

        return QueryValidationInput(
            source=ValidationSource.ASSEMBLY,
            passed=passed,
            lowest_pass_error=lowest_pass_error,
            lowest_error=lowest_error,
            validating_seq=validating_seq,
            segment_results=best_segment_validation_results,
            best_strand_match=best_strand_match,
        )
