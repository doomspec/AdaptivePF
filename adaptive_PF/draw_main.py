from mizore.method.evolution._te_drawer import *
from mizore.draw._standard_plt import get_standard_plt, init_subplots
from draw_compare import draw_compare
from draw_krylov import draw_krylov
from draw_average import draw_average
from mizore.utilities.path_tools import get_subdirectory_paths
import numpy as np


if __name__ == '__main__':
    root_path="/home/mh-group/mizore_results/adaptive_evolution"
    status_path = "/home/mh-group/mizore_results/adaptive_evolution/h4/H4_0.05_01-07-09h24m00s"

    plt = get_standard_plt()

    fig, axs = init_subplots(plt, ncols=3, nrows=2, fig_width=10, wid_high_ratio=0.57, more_options={"auto_title": True})
    draw_compare([axs[0][0],axs[1][0]])
    draw_krylov([axs[0][1]])
    draw_run_status_with_trotter_gate(axs[1][1],status_path)
    draw_average([axs[0][2],axs[1][2]])

    plt.subplots_adjust(hspace=0.2,wspace=0.3)

    plt.savefig(root_path + '/numerics_new.pdf', bbox_inches='tight')
