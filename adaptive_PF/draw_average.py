import json, os
from mizore.utilities.path_tools import get_subdirectory_paths
from mizore.method.evolution._te_drawer import get_fidelity_x_y_data
def folder_name2title(folder_name):
    cutoff=float(folder_name.split("_")[1])
    return "$\Delta_{\mathrm{cut}}="+str(cutoff)+"$"

def trotter_name2title(trotter_name):
    n_delta_t = 2
    divider, divided = trotter_name[7:].split("-")
    if n_delta_t!=1:
        print("Notice: n_delta_t!=1")
    print("Please remove *5 in the next line !!!! draw_compare trotter_name2title")
    return "Trotter {}".format(int(float(divided)/float(divider)*n_delta_t)*5)

def get_max_min_avg_dict(input_lists):
    """
    :return: {"max":[...],"min":[...],"avg":[...]}
    """
    n_point=len(input_lists[0])
    max_list = [None] * n_point
    min_list = [None] * n_point
    avg_list = [None] * n_point
    for i in range(n_point):
        small_list=[]
        for input_list in input_lists:
            small_list.append(input_list[i])
        max_list[i]=max(small_list)
        min_list[i]=min(small_list)
        avg_list[i]=sum(small_list)/len(small_list)

    max_min_avg_list_dict={"max":max_list,"min":min_list,"avg":avg_list}

    return max_min_avg_list_dict

def draw_max_min_avg_dict_and_save(ax,max_min_avg_dict,delta_t):
    #print(max_min_avg_list_dict)
    ax.grid()
    ax.set_xlabel("Time")
    ax.set_ylabel("Infidelity")
    time_list = [delta_t * i for i in range(len(max_min_avg_dict["main"]["avg"]))]
    n_step = len(time_list)
    for key in max_min_avg_dict.keys():

        if key == "main":
            label = "Adaptive PF"
            zorder = 10
        else:
            label = trotter_name2title(key)
            zorder=0

        local_time_list,avg_list=get_fidelity_x_y_data(max_min_avg_dict[key]["avg"],delta_t)
        ax.plot(local_time_list, avg_list, '-',
                label=label,zorder=zorder,alpha=0.9)
        min_list=[num for num in max_min_avg_dict[key]["min"] if num>=0]
        max_list=[num for num in max_min_avg_dict[key]["max"] if num>=0]
        ax.fill_between(local_time_list,min_list,max_list, alpha=0.2)
    ax.legend()
    format_exponent(ax)


def get_delta_t_and_max_min_avg_list_dict_from_paths(paths,key2ignore=tuple()):

    with open(paths[0] + "/para_dict.json", "r") as f:
        info_dict = json.load(f)

    delta_t = info_dict["delta_t"]

    infidelity_dict_list=[]
    for path in paths:
        with open(path + "/benchmark/infidelity.json", "r") as f:
            infidelity_dict = json.load(f)
            infidelity_dict_list.append(infidelity_dict)
            #print(path,sum(infidelity_dict["main"]))

    infidelity_list_dict={key:[] for key in infidelity_dict.keys()}
    for i_dict in infidelity_dict_list:
        for key,i_list in i_dict.items():
            infidelity_list_dict[key].append(i_list)

    max_min_avg_list_dict = {}
    for key,i_lists in infidelity_list_dict.items():
        if key not in key2ignore:
            max_min_avg_list_dict[key]=get_max_min_avg_dict(i_lists)

    return delta_t, max_min_avg_list_dict

from mizore.method.evolution._te_drawer import *
from mizore.draw._standard_plt import get_standard_plt, init_subplots


def draw_gate_compare_same_color(ax, paths, gate_name="CNOT"):
    gate_num_dict = {}
    gate_time_dict = {}
    for project_path in paths:
        gate_num_list = []
        gate_time_list = []
        with open("{}/run_status_info.json".format(project_path), "r") as f:
            log_dict = json.load(f)

        gate_count_list = log_dict["n_block_change"]
        for i_gate_count in range(len(log_dict["n_block_change"])):
            gate_num_list.append(gate_count_list[i_gate_count][2][gate_name])
            gate_time_list.append(gate_count_list[i_gate_count][0])
            if i_gate_count != len(log_dict["n_block_change"]) - 1:
                gate_num_list.append(gate_count_list[i_gate_count][2][gate_name])
                gate_time_list.append(gate_count_list[i_gate_count + 1][0])

        with open("{}/para_dict.json".format(project_path), "r") as f:
            para_dict = json.load(f)
        gate_num_list.append(gate_num_list[-1])
        gate_time_list.append(para_dict["n_circuit"] * para_dict["delta_t"])
        folder_name = project_path.split("/")[-1]
        gate_num_dict[folder_name] = gate_num_list
        gate_time_dict[folder_name] = gate_time_list

    for key in gate_time_dict.keys():
        ax.plot(gate_time_dict[key], gate_num_dict[key], "-", label=key, color="steelblue", alpha=0.6)
    format_exponent(ax)
    ax.legend()
    ax.grid()
    ax.set_xlabel("Time")
    ax.set_ylabel(gate_name + " gate count")

folder_path="/home/mh-group/mizore_results/adaptive_evolution/random_ising"

def draw_average(axs):
    paths = get_subdirectory_paths(folder_path)
    print("N_projects:",len(paths))
    key2ignore=["trotter1-1"]
    delta_t, max_min_avg_list_dict=get_delta_t_and_max_min_avg_list_dict_from_paths(paths,key2ignore=key2ignore)
    draw_max_min_avg_dict_and_save(axs[0],max_min_avg_list_dict,delta_t)
    draw_gate_compare_same_color(axs[1],paths)
    axs[1].get_legend().remove()

if __name__ == '__main__0':
    from benchmark import benchmark
    paths = get_subdirectory_paths(folder_path)
    for name in paths:
        benchmark(name,trotter=[{"n_delta_t": 2, "n_trotter_step": 3}],exact_state_generated=True)

if __name__ == '__main__':
    plt = get_standard_plt()
    fig, axs = init_subplots(plt, ncols=2, nrows=1)
    draw_averge(axs)
    plt.savefig(folder_path + '/compare.png', bbox_inches='tight')