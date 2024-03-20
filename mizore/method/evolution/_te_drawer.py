import pickle, json
from ._utilities import get_first_trotter_gate_count

from mizore.draw.utilities import format_exponent

def draw_run_status(ax, path, gate_name="CNOT"):

    with open(path + "/hamiltonian.pickle", "rb") as f:
        hamiltonian = pickle.load(f)
    with open(path + "/run_status_info.json", "r") as f:
        log_dict = json.load(f)
    n_block_change = log_dict["n_block_change"]
    n_gate_list = []
    n_gate_time_list = []
    for item in n_block_change:
        n_gate_time_list.append(item[0])
        n_gate_list.append(item[2][gate_name])
    quality_list = log_dict["quality_list"]
    evolution_time_list = log_dict["evolution_time_list"]
    first_trotter_gate_use = get_first_trotter_gate_count(hamiltonian,gate_name=gate_name)
    end_time = evolution_time_list[-1]

    n_gate_list = []
    n_gate_time_list = []
    n_block_change.append((end_time, None, None))
    for i in range(len(n_block_change) - 1):
        n_gate_time_list.append(n_block_change[i][0])
        n_gate_list.append(n_block_change[i][2][gate_name])
        n_gate_time_list.append(n_block_change[i + 1][0])
        n_gate_list.append(n_block_change[i][2][gate_name])

    lns1 = ax.plot(evolution_time_list, quality_list, '-',
                   color="dodgerblue", label="Distance")
    ax2 = ax.twinx()
    lns2 = ax2.plot(n_gate_time_list, n_gate_list, '-',
                    color="orange", label=gate_name + " count: Adaptive")
    #lns3 = ax2.plot([0, end_time], [first_trotter_gate_use, first_trotter_gate_use],
     #               '-', color="red", label=gate_name + " count: 1st Trotter ")
    ax.grid()
    ax.set_xlabel("Time")
    ax.set_ylabel("Distance")
    ax2.set_ylabel("Number of gates")
    lns = lns2 + lns1 #+ lns3
    labs = [l.get_label() for l in lns]
    ax.legend(lns, labs, loc='upper right', ncol=2)


    gate_used_lim = n_gate_list[-1] # max(first_trotter_gate_use, n_gate_list[-1])
    ax2.set_ylim(0, gate_used_lim * 1.3)
    ax.set_ylim(0, max(quality_list) * 1.4)
    format_exponent(ax)

def draw_run_status_with_trotter_gate(ax, path, gate_name="CNOT"):

    with open(path + "/hamiltonian.pickle", "rb") as f:
        hamiltonian = pickle.load(f)
    with open(path + "/run_status_info.json", "r") as f:
        log_dict = json.load(f)
    n_block_change = log_dict["n_block_change"]
    n_gate_list = []
    n_gate_time_list = []
    for item in n_block_change:
        n_gate_time_list.append(item[0])
        n_gate_list.append(item[2][gate_name])
    quality_list = log_dict["quality_list"]
    evolution_time_list = log_dict["evolution_time_list"]
    first_trotter_gate_use = get_first_trotter_gate_count(hamiltonian,gate_name=gate_name)
    end_time = evolution_time_list[-1]

    n_gate_list = []
    n_gate_time_list = []
    n_block_change.append((end_time, None, None))
    for i in range(len(n_block_change) - 1):
        n_gate_time_list.append(n_block_change[i][0])
        n_gate_list.append(n_block_change[i][2][gate_name])
        n_gate_time_list.append(n_block_change[i + 1][0])
        n_gate_list.append(n_block_change[i][2][gate_name])

    lns1 = ax.plot(evolution_time_list, quality_list, '-',
                   color="dodgerblue", label="Error $\Delta$")
    ax2 = ax.twinx()
    lns2 = ax2.plot(n_gate_time_list, n_gate_list, '-',
                    color="orange", label="# "+gate_name + ": Adaptive PF")
    lns3 = ax2.plot([0, end_time], [first_trotter_gate_use, first_trotter_gate_use],
                    '-', color="red", label="# "+gate_name + ": A Trotter step")
    ax.grid()
    ax.set_xlabel("Time")
    ax.set_ylabel("Error $\Delta$")
    ax2.set_ylabel(gate_name+" gate count")
    ax2.ticklabel_format(axis="y", style='sci', scilimits=(-2, 2),useMathText=True)
    lns = lns3 + lns2 + lns1
    labs = [l.get_label() for l in lns]
    ax.legend(lns, labs, loc="upper left", ncol=2)


    gate_used_lim = max(first_trotter_gate_use, n_gate_list[-1])
    ax2.set_ylim(0, gate_used_lim * 1.3)
    _ , y_max=ax.get_ylim()
    ax.set_ylim(0, y_max * 1.45)
    format_exponent(ax)


def draw_benchmark(ax,path, key2ignore=(None,)):

    with open(path+"/benchmark/infidelity.json","r") as f:
        infidelity_dict=json.load(f)
    with open(path+"/para_dict.json","r") as f:
        delta_t=json.load(f)["delta_t"]

    max_infidelity=0

    for label in infidelity_dict.keys():
        if label in key2ignore:
            continue
        infidelity_list = infidelity_dict[label]
        #x_data = [step * delta_t for step in range(len(infidelity_list))]
        x_data, y_data = get_fidelity_x_y_data(infidelity_list, delta_t)
        ax.plot(x_data, y_data, '-o',label=label, alpha=0.8)
        local_max_infidelity = max(infidelity_list)
        if local_max_infidelity > max_infidelity:
            max_infidelity = local_max_infidelity

    ax.grid()
    ax.set_xlabel("Time")
    ax.set_ylabel("1-Fidelity")
    ax.legend(loc='upper right', ncol=3)
    max_diff = max_infidelity
    ax.set_ylim(-max_diff*0.05, max_diff*1.2)
    format_exponent(ax)

def get_fidelity_x_y_data(infidelity_list, delta_t):
    filtered_list=[]
    t_list=[]
    for step in range(len(infidelity_list)):
        if infidelity_list[step]>=0:
            filtered_list.append(infidelity_list[step])
            t_list.append(step*delta_t)
    return t_list,filtered_list