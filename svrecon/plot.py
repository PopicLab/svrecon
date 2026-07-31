import logging
from pathlib import Path
from typing import Optional, Sequence, TYPE_CHECKING

import matplotlib.pyplot as plt
import wotplot as wp
from matplotlib.ticker import FuncFormatter, MaxNLocator

if TYPE_CHECKING:
    from svrecon.scoring import QueryInfo

logger = logging.getLogger(__name__)

# Plotting Hyperparameters
k = 15
TICK_COUNT = 8  # target number of ticks per axis, regardless of sequence length

# Functions
def plot_dot_plot(s1: str, s2: str, title: str, output_path: str,
                  s1_name: str = 'reconstructed', s2_name: str = 'validating',
                  segment_boundaries: Optional[Sequence[int]] = None,
                  x_offset: int = 0, y_offset: int = 0) -> None:
    """x_offset/y_offset shift the tick LABELS (not the data) so an axis whose sequence
    starts at a known genomic position can display real coordinates instead of raw
    0-based indices into s1/s2 -- e.g. y_offset=ref_start when s2 is a reference slice."""
    matrix = wp.DotPlotMatrix(s1.upper(), s2.upper(), k)

    # s1_name labels the x-axis, s2_name the y-axis (wotplot's own convention).
    fig, ax = wp.viz_imshow(matrix, title=title, s1_name=s1_name, s2_name=s2_name)

    # Nice, evenly-spaced ticks that scale to the sequence length instead of wotplot's
    # default density, so a 200 bp window and a 50 kb window are equally readable.
    ax.xaxis.set_major_locator(MaxNLocator(nbins=TICK_COUNT, integer=True))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=TICK_COUNT, integer=True))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _pos: f'{int(x) + x_offset:,}'))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _pos: f'{int(y) + y_offset:,}'))

    if segment_boundaries:
        # Boundaries are in s1 (x-axis) coordinates -- mark where each reconstructed
        # piece (junction/flank) begins or ends.
        for pos in segment_boundaries:
            ax.axvline(x=pos, color='red', linestyle='--', linewidth=0.6, alpha=0.6)

    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_query_dot_plots(query_info: 'QueryInfo', svid: str, sv_type: str, output_dir: str, part: int = 0) -> None:
    """Writes two dot-plot PNGs for one reconstructed subsequence of an SV (`part` is that
    subsequence's index among the SV's possibly-compound set of reconstructed pieces):
    1. reconstructed vs. validating sequence -- the real read/assembly/edlib span that confirmed it.
    2. reconstructed vs. reference sequence -- the unmodified reference window it was built from.
    Segment boundaries (query_info.query.segments) are drawn as vertical lines on both, marking
    where each reconstructed piece (junction, flank) begins/ends. The title carries the SV's
    genomic locus (chrom:ref_start-ref_end) since the reconstructed sequence's own coordinates
    don't map 1:1 to the genome once pieces have been rearranged."""
    recon_seq = query_info.query.sequence
    boundaries = sorted({pos for start, end in query_info.query.segments for pos in (start, end)})
    source = query_info.source.value if query_info.source else 'none'

    chrom = query_info.query.chrom
    ref_start = query_info.query.ref_start
    ref_end = query_info.query.ref_end
    locus = f'{chrom}:{ref_start:,}-{ref_end:,}'
    title = f'{sv_type} {svid}_{part} ({locus})'

    if query_info.validating_seq is not None:
        plot_path = Path(output_dir) / f'{sv_type}_{svid}_{part}_validating_{source}.png'
        plot_dot_plot(recon_seq, query_info.validating_seq, title, str(plot_path),
                     s1_name='reconstructed', s2_name=f'validating ({source})',
                     segment_boundaries=boundaries)
    else:
        logger.warning(f'Skipping reconstructed-vs-validating dot plot for {svid}_{part} ({sv_type}): no validating sequence.')

    if query_info.ref_sequence is not None:
        plot_path = Path(output_dir) / f'{sv_type}_{svid}_{part}_reference.png'
        # The reference slice starts at ref_start, so its own position 0 IS genomic
        # position ref_start -- offset the y-axis ticks to show real coordinates.
        plot_dot_plot(recon_seq, query_info.ref_sequence, title, str(plot_path),
                     s1_name='reconstructed', s2_name=f'reference',
                     segment_boundaries=boundaries, y_offset=ref_start)
    else:
        logger.warning(f'Skipping reconstructed-vs-reference dot plot for {svid}_{part} ({sv_type}): no reference sequence.')
