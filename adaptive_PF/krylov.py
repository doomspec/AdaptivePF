from mizore.method.subspace._eigensolver import SubspaceEigensolver
from mizore.method.evolution._te_benchmarker import *
from mizore.classical.numpy_core import get_kth_excited_state
from openfermion.utils import count_qubits
import pickle
from mizore.utilities.path_tools import get_subdirectory_paths



def run_krylov(save_path, n_delta_t=1, key2run=None):

    if key2run is None:
        key2run = ["main"]

    with open(save_path + "/hamiltonian.pickle", "rb") as f:
        hamiltonian = pickle.load(f)

    state_dict = get_project_quantum_state_dict(save_path, n_delta_t=n_delta_t,key2run=key2run)

    solver = SubspaceEigensolver()
    solver.load_hamiltonian_ops(hamiltonian)

    """
    n_qubit=count_qubits(hamiltonian)
    fci_energy, eigvec = get_kth_excited_state(0, n_qubit, hamiltonian)
    """


    ground_energy_dict={}

    for key, value in state_dict.items():
        solver.load_quantum_states(value)
        ground_energy = solver.solve()
        ground_energy_dict[key]=ground_energy

    """
    solver.load_quantum_states([value[0]])
    ground_energy = solver.solve()
    ground_energy_dict["hf"] = ground_energy
    """

    return ground_energy_dict

def get_fci_by_path(save_path):
    with open(save_path + "/hamiltonian.pickle", "rb") as f:
        hamiltonian = pickle.load(f)
    n_qubit = count_qubits(hamiltonian)
    fci_energy, eigvec = get_kth_excited_state(0, n_qubit, hamiltonian)
    return fci_energy

def get_n_gate_by_path(save_path):
    with open(save_path + "/run_status_info.json", "r") as f:
        run_status_dict = json.load(f)
    return run_status_dict["n_block_change"][-1][-1]["CNOT"]

from benchmark import benchmark
from draw_single import draw_single

root_path = "/home/mh-group/mizore_results/adaptive_evolution/h4"
paths = get_subdirectory_paths(root_path)

if __name__ == '__main__0':
    for name in paths:
        benchmark(name,trotter=[{"n_delta_t": 2, "n_trotter_step": 1}],exact_state_generated=True)

if __name__ == '__main__':
    n_delta_t = 2
    res_dict={"run":{},"benchmark":{},"setting":{}, "n_gate":{}}
    for name in paths:
        d_cut=float(name.split("/")[-1].split("_")[1])
        res_dict["run"][d_cut]=run_krylov(name, n_delta_t=n_delta_t, key2run=["main"])["main"]
        res_dict["n_gate"][d_cut]=get_n_gate_by_path(name)

    res_dict["benchmark"].update(run_krylov(paths[0], n_delta_t=n_delta_t, key2run=["trotter1-1","trotter2-1","exact"]))
    res_dict["benchmark"]["fci"]=get_fci_by_path(paths[0])
    res_dict["setting"]["n_delta_t"]=n_delta_t

    with open(root_path+"/krylov.json","w") as f:
        json.dump(res_dict,f)
    print(res_dict)
