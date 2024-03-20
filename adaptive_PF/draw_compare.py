from mizore.method.evolution._te_drawer import *
from mizore.draw._standard_plt import get_standard_plt, init_subplots
from mizore.utilities.path_tools import get_subdirectory_paths

def folder_name2title(folder_name):
    cutoff=float(folder_name.split("_")[1])
    return "$\Delta_{\mathrm{cut}}="+str(cutoff)+"$"

def trotter_name2title(trotter_name):
    divider, divided = trotter_name[7:].split("-")
    if n_delta_t!=1:
        print("Notice: n_delta_t!=1")
    print("Please remove *15 in the next line !!!! draw_compare trotter_name2title")
    return "Trotter {}".format(int(float(divided)/float(divider)*n_delta_t)*15)

def draw_gate_compare(ax, paths, gate_name="CNOT"):

    gate_num_dict = {}
    gate_time_dict = {}
    for project_path in paths:
        gate_num_list = []
        gate_time_list = []
        with open("{}/run_status_info.json".format(project_path), "r") as f:
            log_dict = json.load(f)

        gate_count_list=log_dict["n_block_change"]
        for i_gate_count in range(len(log_dict["n_block_change"])):
            gate_num_list.append(gate_count_list[i_gate_count][2][gate_name])
            gate_time_list.append(gate_count_list[i_gate_count][0])
            if i_gate_count!=len(log_dict["n_block_change"])-1:
                gate_num_list.append(gate_count_list[i_gate_count][2][gate_name])
                gate_time_list.append(gate_count_list[i_gate_count+1][0])

        with open("{}/para_dict.json".format(project_path), "r") as f:
            para_dict = json.load(f)
        gate_num_list.append(gate_num_list[-1])
        gate_time_list.append(para_dict["n_circuit"] * para_dict["delta_t"])
        folder_name=project_path.split("/")[-1]
        gate_num_dict[folder_name] = gate_num_list
        gate_time_dict[folder_name] = gate_time_list

    for key in gate_time_dict.keys():
        ax.plot(gate_time_dict[key], gate_num_dict[key], "-", label=folder_name2title(key))

    ax.legend(ncol=3,loc="upper left")
    ax.grid()
    ax.set_xlabel("Time")
    ax.set_ylabel(gate_name + " gate count")
    ax.ticklabel_format(axis="y", style='sci', scilimits=(-2, 2), useMathText=True)
    ymin,ymax=ax.get_ylim()
    ax.set_ylim(ymin,ymax*1.3)


def draw_fidelity_compare(ax, paths, benchmark2draw=tuple()):
    infidelity_dict = {}
    for name in paths:
        with open("{}/benchmark/infidelity.json".format(name), "r") as f:
            saved_infidelity_dict = json.load(f)
        infidelity_dict[name] = saved_infidelity_dict["main"]


    with open("{}/para_dict.json".format(paths[0]), "r") as f:
        delta_t = json.load(f)["delta_t"]

    max_infidelity=-1

    for key in infidelity_dict.keys():
        x_data, y_data = get_fidelity_x_y_data(infidelity_dict[key], delta_t)
        folder_name = key.split("/")[-1]
        ax.plot(x_data, y_data, "-", label=folder_name2title(folder_name))
        max_infidelity=max(max_infidelity,max(y_data))

    legend_lns=None
    legend_labs = []
    for key in benchmark2draw:
        #print(saved_infidelity_dict)
        infidelity4benchmark = saved_infidelity_dict[key]
        x_data, y_data = get_fidelity_x_y_data(infidelity4benchmark, delta_t)
        lns=ax.plot(x_data, y_data, "-", linewidth=5,alpha=0.5)
        if legend_lns is None:
            legend_lns=lns
        else:
            legend_lns=legend_lns+lns
        legend_labs.append(trotter_name2title(key))
    if legend_labs:
        ax.legend(legend_lns,legend_labs,ncol=4,loc = 'upper left',columnspacing=0.8,fontsize=15)
    ax.grid()
    ax.set_ylim(0,max_infidelity*1.3)
    ax.set_xlabel("Time")
    ax.ticklabel_format(axis="y", style='sci', scilimits=(-2, 2), useMathText=True)
    ax.set_ylabel("Infidelity")

def order_path_by_dcut(paths):
    d_cut_list=[]
    for name in paths:
        d_cut_list.append(float(name.split("/")[-1].split("_")[1]))
    from numpy import argsort
    rank=argsort(d_cut_list)
    paths=[paths[i] for i in rank]
    return  paths

#root_path = "/home/mh-group/mizore_results/adaptive_evolution/h4"
root_path = "/home/mh-group/mizore_results/adaptive_evolution/h2o"
n_delta_t=2
paths = get_subdirectory_paths(root_path)
paths = order_path_by_dcut(paths)


def draw_compare_bak(axs):
    draw_gate_compare(axs[0], paths)
    draw_fidelity_compare(axs[1], paths)

def draw_compare(axs):
    draw_gate_compare(axs[0], paths)
    draw_fidelity_compare(axs[1], paths, benchmark2draw=["trotter1-1","trotter2-3","trotter1-2","trotter1-3"])

from benchmark import benchmark

if __name__ == '__main__0':
    for name in paths:
        benchmark(name,trotter=[{"n_delta_t": 3, "n_trotter_step": 4}],exact_state_generated=True)

if __name__ == '__main__':
    plt = get_standard_plt()
    fig, axs = init_subplots(plt, ncols=2, nrows=1)
    draw_compare(axs)
    plt.savefig(root_path + '/compare.png', bbox_inches='tight')