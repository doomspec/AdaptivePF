from mizore.method.evolution._te_drawer import *
from mizore.draw._standard_plt import get_standard_plt, init_subplots
from mizore.draw import format_exponent
from mizore.utilities.path_tools import get_subdirectory_paths
#from draw_utils import *
from numpy import argsort

def trotter_name2title_with_n_delta_t(trotter_name,n_delta_t):
    if trotter_name=="exact":
        return "Exact"
    divider, divided = trotter_name[7:].split("-")
    print("Please remove *15 in the next line !!!! draw_compare trotter_name2title")
    return "Trotter {}".format(int(float(divided)/float(divider)*n_delta_t)*15)

def draw_krylov_compare(ax, json_path, benchmark2draw=tuple()):

    with open(json_path, "r") as f:
        krylov_dict = json.load(f)

    fci_energy = krylov_dict["benchmark"]["fci"]
    n_delta_t=krylov_dict["setting"]["n_delta_t"]

    x_data=[]
    energy_data=[]
    n_gate_data=[]
    for key, value in krylov_dict["run"].items():
        x_data.append(float(key))
        energy_data.append(value)
        n_gate_data.append(krylov_dict["n_gate"][key])

    order_index = argsort(x_data)
    x_data = [x_data[i] for i in order_index]
    energy_data = [energy_data[i]-fci_energy for i in order_index]
    n_gate_data = [n_gate_data[i] for i in order_index]
    print(energy_data)
    ax2=ax.twinx()

    labels=["Adaptive PF: Energy"]
    lns=ax.plot(x_data, energy_data, "-o")
    lns += ax2.plot(x_data, n_gate_data, "-s", color="orange")
    labels.append("Adaptive PF: # CNOT")

    for name in benchmark2draw:
        value=krylov_dict["benchmark"][name]-fci_energy
        lns+=ax.plot([x_data[0], x_data[-1]], [value, value], "--")
        labels.append(trotter_name2title_with_n_delta_t(name,n_delta_t))
        pass

    chemical_accuracy = 1e-3
    ax.plot([x_data[0], x_data[-1]], [chemical_accuracy, chemical_accuracy], "-", color="red",zorder=-1)
    #labels.append("Chemical accuracy")



    ax.invert_xaxis()
    ax.grid()
    ax.set_xlabel("Error cutoff $\Delta_{\mathrm{cut}}$")
    ax.set_ylabel("Energy error [Hartree]")
    ax2.set_ylabel("CNOT gate count")
    ax2.ticklabel_format(axis="y", style='sci', scilimits=(-2, 2), useMathText=True)
    format_exponent(ax)
    ax.set_ylim(ax.get_ylim()[0],max(energy_data) * 1.3)
    ax2.set_ylim(min(n_gate_data) * 0.95,max(n_gate_data)*1.2)
    ax.legend(lns,labels,ncol=3,loc="upper left")

root_path = "/home/mh-group/mizore_results/adaptive_evolution/h4"
paths = get_subdirectory_paths(root_path)

def draw_krylov(axs):
    draw_krylov_compare(axs[0], root_path+"/krylov.json",benchmark2draw=["trotter2-1","trotter1-1","exact"])

if __name__ == '__main__':
    plt = get_standard_plt()
    fig, axs = init_subplots(plt, ncols=1, nrows=1)
    draw_krylov(axs)
    plt.savefig(root_path + '/krylov.png', bbox_inches='tight')