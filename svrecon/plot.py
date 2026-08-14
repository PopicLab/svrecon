import logging
from pathlib import Path
from typing import Optional, Sequence, TYPE_CHECKING

import matplotlib.pyplot as plt
import wotplot as wp
from matplotlib.ticker import FuncFormatter, MaxNLocator

if TYPE_CHECKING:
    from svrecon.scoring import QueryInfo

logger = logging.getLogger(__name__)

K = 15  # wotplot k-mer size
TICK_COUNT = 8  # target ticks per axis, regardless of sequence length


def plot_dot_plot(s1: str, s2: str, title: str, output_path: str,
                  s1_name: str = 'reconstructed', s2_name: str = 'validating',
                  segment_boundaries: Optional[Sequence[int]] = None,
                  ref_segment_boundaries: Optional[Sequence[int]] = None,
                  x_offset: int = 0, y_offset: int = 0, aspect: str = 'auto') -> None:
    matrix = wp.DotPlotMatrix(s1.upper(), s2.upper(), K)
    num_rows = matrix.mat.shape[0]
    fig, ax = wp.viz_spy(matrix, markersize=1.0, aspect=aspect, title=title, s1_name=s1_name, s2_name=s2_name)
    ax.xaxis.tick_bottom()  # spy() defaults to top-side ticks; move them to match the x-label below

    # Ticks scale with sequence length instead of wotplot's default density.
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


def plot_query_dot_plots(query_info: 'QueryInfo', svid: str, sv_type: str, output_dir: str, part: int = 0,
                         aspect: str = 'auto') -> None:
    """Writes reconstructed-vs-validating and reconstructed-vs-reference dot plots for one
    subsequence directly into output_dir, which must already exist -- the caller (score_all)
    owns the directory layout (organizing by validation outcome, SV type, and SVID) and its
    creation, since that's a scoring concern, not a plotting one. Ref-coordinate boundaries
    are drawn only on the latter, since a piece's ref position can diverge from its alt
    position once moved, inverted, or pasted."""
    query = query_info.query
    recon_seq = query.sequence
    alt_boundaries = sorted({pos for start, end in query.recon_segments for pos in (start, end)})
    source = query_info.source.value if query_info.source else 'none'
    locus = f'{query.chrom}:{query.ref_start:,}-{query.ref_end:,}'
    title = f'{sv_type} {svid}_{part} ({locus})'
    out_dir = Path(output_dir)

    def write(seq, plot_suffix, s2_name, seq_label, **extra):
        if seq is None:
            logger.info(f'Skipping reconstructed-vs-{seq_label} dot plot for {svid}_{part} ({sv_type}): '
                       f'no {seq_label} sequence.')
            return
        plot_path = out_dir / f'{part}_{plot_suffix}.png'
        plot_dot_plot(recon_seq, seq, title, str(plot_path), s1_name='reconstructed', s2_name=s2_name,
                     segment_boundaries=alt_boundaries, aspect=aspect, **extra)

    write(query_info.validating_seq, f'validating_{source}', f'validating ({source})', 'validating')

    ref_boundaries = sorted({pos - query.ref_start for start, end in query.ref_segments for pos in (start, end)})
    write(query_info.ref_sequence, 'reference', 'reference', 'reference',
         ref_segment_boundaries=ref_boundaries, y_offset=query.ref_start)
