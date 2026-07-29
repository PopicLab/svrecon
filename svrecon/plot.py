import matplotlib.pyplot as plt
import wotplot as wp

# Plotting Hyperparameters
k = 15

# Functions
def plot_dot_plot(s1: str, s2: str, title: str, output_path: str,
                  s1_name: str = 'reconstructed', s2_name: str = 'validating') -> None:
    matrix = wp.DotPlotMatrix(s1.upper(), s2.upper(), k)

    # s1_name labels the x-axis, s2_name the y-axis (wotplot's own convention).
    fig, ax = wp.viz_imshow(matrix, title=title, s1_name=s1_name, s2_name=s2_name)
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
