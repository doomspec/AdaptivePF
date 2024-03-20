from mizore.method.evolution._te_drawer import *
from mizore.draw._standard_plt import get_standard_plt, init_subplots


def draw_single(project_path):
    plt = get_standard_plt()
    fig, axs = init_subplots(plt, ncols=2, nrows=1)
    draw_benchmark(axs[1], project_path, key2ignore=("init",))
    draw_run_status(axs[0], project_path)
    plt.savefig(project_path + '/run_status.png', bbox_inches='tight')

if __name__ == '__main__':
    draw_single("/home/mh-group/mizore_results/adaptive_evolution/H2O_0.2_01-07-06h09m38s")
