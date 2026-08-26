import logging
import random
import re
from pathlib import Path
from typing import Optional, Sequence, TYPE_CHECKING

import matplotlib.pyplot as plt
import wotplot as wp
from matplotlib.ticker import FuncFormatter, MaxNLocator

if TYPE_CHECKING:
    from svrecon.scoring import SVValidationResult

logger = logging.getLogger(__name__)

K = 15  # wotplot k-mer size
TICK_COUNT = 8  # target ticks per axis, regardless of sequence length
AXIS_FONTSIZE = 18 
PLOTTABLE_BASES = 'ACGT'  # wotplot rejects N and every other IUPAC code
UNPLOTTABLE_BASE = re.compile(f'[^{PLOTTABLE_BASES}]')


def plot_dot_plot(s1: str, s2: str, title: str, output_path: str,
                  s1_name: str = 'reconstructed', s2_name: str = 'validating',
                  segment_boundaries: Optional[Sequence[int]] = None,
                  ref_segment_boundaries: Optional[Sequence[int]] = None,
                  x_offset: int = 0, y_offset: int = 0, aspect: str = 'auto',
                  axis_length: bool = False) -> None:
    """Renders a k-mer dot plot of s1 (x) vs s2 (y) to output_path, with optional dashed
    boundary lines per axis and x/y tick offsets for genomic coordinates."""
    matrix = wp.DotPlotMatrix(s1.upper(), s2.upper(), K)

    num_rows = matrix.mat.shape[0]
    fig, ax = wp.viz_spy(matrix, markersize=1.0, aspect=aspect, title=title, s1_name=s1_name, s2_name=s2_name)     # Ticks scale with sequence length instead of wotplot's default density.
    ax.set_xlabel(f'{s1_name} ({len(s1):,} nt)' if axis_length else s1_name, fontsize=AXIS_FONTSIZE)
    ax.set_ylabel(f'{s2_name} ({len(s2):,} nt)' if axis_length else s2_name, fontsize=AXIS_FONTSIZE)
    ax.xaxis.tick_bottom()  # spy() defaults to top-side ticks; move them to match the x-label below
    ax.xaxis.set_major_locator(MaxNLocator(nbins=TICK_COUNT, integer=True))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=TICK_COUNT, integer=True))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _pos: f'{int(x) + x_offset:,}'))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _pos: f'{num_rows - 1 - int(y) + y_offset:,}'))

    # segment_boundaries: s1/x-axis coords. ref_segment_boundaries: s2/y-axis, row-indexed like the ticks.
    for pos in segment_boundaries or ():
        ax.axvline(x=pos, color='gray', linestyle='--', linewidth=0.6, alpha=0.6)
    for pos in ref_segment_boundaries or ():
        ax.axhline(y=num_rows - 1 - pos, color='blue', linestyle='--', linewidth=0.6, alpha=0.6)

    fig.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close(fig)


def plot_sv_validation(sv_validation: 'SVValidationResult', output_dir: str, aspect: str = 'auto',
                       substitute_bases: bool = False, title_svid: bool = False,
                       title_location: bool = False, axis_length: bool = False) -> None:
    """Plots each subsequence of one SV vs the reference and (when found) its matched
    sequence, into an existing output_dir -- the caller owns directory layout and gating."""
    out_dir = Path(output_dir)
    for part, query_validation in enumerate(sv_validation.query_validation_results):
        query = query_validation.query
        alt_boundaries = sorted({pos for start, end in query.recon_segments for pos in (start, end)})
        title = query.svtype
        if title_svid:
            title += f' {query.svid}_{part}'
        if title_location:
            title += f' ({query.chrom}:{query.ref_start:,}-{query.ref_end:,})'

        recon_seq = query.sequence.upper()
        ref_seq = query.ref_sequence.upper()
        matched_seq = (query_validation.best_matched_seq or '').upper()

        # skip plotting if non plottable base pair exists in sequence, unless substituted
        plottable = not UNPLOTTABLE_BASE.search(recon_seq + ref_seq + matched_seq)
        if not plottable and substitute_bases:
            # substitute unknown bases with random ATCG
            rng = random.Random(0)  # seeded so a re-run draws the same bases
            random_base = lambda _match: rng.choice(PLOTTABLE_BASES)
            recon_seq, num_recon_subs = UNPLOTTABLE_BASE.subn(random_base, recon_seq)
            ref_seq, num_ref_subs = UNPLOTTABLE_BASE.subn(random_base, ref_seq)
            matched_seq, num_matched_subs = UNPLOTTABLE_BASE.subn(random_base, matched_seq)
            logger.info(f'Substituting {num_recon_subs + num_ref_subs + num_matched_subs} '
                        f'bases for {query.svid}')
            plottable = True
        if not plottable:
            logger.info(f'Skipping {query.svid}_{part} plots, sequence contains unknown basepair')
            continue

        # reconstructed vs the unmodified reference
        ref_boundaries = sorted({pos - query.ref_start
                                 for start, end in query.ref_segments for pos in (start, end)})
        plot_dot_plot(recon_seq, ref_seq, title, str(out_dir / f'{part}_reference.png'),
                      s2_name='reference', segment_boundaries=alt_boundaries,
                      ref_segment_boundaries=ref_boundaries, y_offset=query.ref_start, aspect=aspect,
                      axis_length=axis_length)

        if matched_seq:
            source = query_validation.source.value
            # reconstructed vs the matched sequence, when a scorer found one (a pass or a match)
            plot_dot_plot(recon_seq, matched_seq, title,
                        str(out_dir / f'{part}_matched_{source}.png'),
                        s2_name=source, segment_boundaries=alt_boundaries, aspect=aspect,
                        axis_length=axis_length)

            # matched vs the reference:
            plot_dot_plot(matched_seq, ref_seq, title,
                        str(out_dir / f'{part}_matched_{source}_reference.png'),
                        s1_name=source, s2_name='reference',
                        ref_segment_boundaries=ref_boundaries, y_offset=query.ref_start, aspect=aspect,
                        axis_length=axis_length)