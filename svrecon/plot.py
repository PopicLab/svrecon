import logging
from pathlib import Path
from typing import Optional, Sequence, TYPE_CHECKING

import matplotlib.pyplot as plt
import wotplot as wp
from matplotlib.ticker import FuncFormatter, MaxNLocator

if TYPE_CHECKING:
    from svrecon.scoring import SVValidation

logger = logging.getLogger(__name__)

K = 15  # wotplot k-mer size
TICK_COUNT = 8  # target ticks per axis, regardless of sequence length
PLOTTABLE_BASES = {'A', 'C', 'G', 'T'}  # wotplot rejects N and every other IUPAC code


def plot_dot_plot(s1: str, s2: str, title: str, output_path: str,
                  s1_name: str = 'reconstructed', s2_name: str = 'validating',
                  segment_boundaries: Optional[Sequence[int]] = None,
                  ref_segment_boundaries: Optional[Sequence[int]] = None,
                  x_offset: int = 0, y_offset: int = 0, aspect: str = 'auto') -> None:
    """Renders a k-mer dot plot of s1 (x) vs s2 (y) to output_path, with optional dashed
    boundary lines per axis and x/y tick offsets for genomic coordinates."""
    matrix = wp.DotPlotMatrix(s1.upper(), s2.upper(), K)

    num_rows = matrix.mat.shape[0]
    fig, ax = wp.viz_spy(matrix, markersize=1.0, aspect=aspect, title=title, s1_name=s1_name, s2_name=s2_name)     # Ticks scale with sequence length instead of wotplot's default density.
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


def plot_sv_validation(sv_validation: 'SVValidation', output_dir: str, aspect: str = 'auto') -> None:
    """Plots each subsequence of one SV vs the reference and (when found) its validating
    sequence, into an existing output_dir -- the caller owns directory layout and gating."""
    out_dir = Path(output_dir)
    for part, query_validation in enumerate(sv_validation.query_validations):
        query = query_validation.query
        alt_boundaries = sorted({pos for start, end in query.recon_segments for pos in (start, end)})
        locus = f'{query.chrom}:{query.ref_start:,}-{query.ref_end:,}'
        title = f'{query.svtype} {query.svid}_{part} ({locus})'

        # skip plotting if non plottable base pair exists in sequence
        recon_plottable = all(char in PLOTTABLE_BASES for char in query.sequence.upper())

        # reconstructed vs the unmodified reference
        ref_boundaries = sorted({pos - query.ref_start
                                 for start, end in query.ref_segments for pos in (start, end)})
        if recon_plottable and all(char in PLOTTABLE_BASES for char in query.ref_sequence.upper()):
            plot_dot_plot(query.sequence, query.ref_sequence, title, str(out_dir / f'{part}_reference.png'),
                          s2_name='reference', segment_boundaries=alt_boundaries,
                          ref_segment_boundaries=ref_boundaries, y_offset=query.ref_start, aspect=aspect)
        else:
            logger.info(f'Skipping {query.svid}_{part} reference plot, sequence contains unknown basepair')

        # reconstructed vs the validating sequence, when a scorer found one
        if query_validation.validating_seq is not None:
            source = query_validation.source.value
            if recon_plottable and all(char in PLOTTABLE_BASES
                                       for char in query_validation.validating_seq.upper()):
                plot_dot_plot(query.sequence, query_validation.validating_seq, title,
                            str(out_dir / f'{part}_validating_{source}.png'),
                            s2_name=f'validating ({source})', segment_boundaries=alt_boundaries, aspect=aspect)
            else:
                logger.info(f'Skipping {query.svid}_{part} validating plot, sequence contains unknown basepair')