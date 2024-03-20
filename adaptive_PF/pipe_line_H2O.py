from mizore.hamiltonian.TestHamiltonian import make_example_12_spin_orbital_H2O, make_example_H2O
from mizore.method.evolution import run_adaptive_product_formula

from benchmark import benchmark
from draw_single import draw_single

if __name__ == '__main__':
    obj = make_example_12_spin_orbital_H2O()
    print(obj.n_qubit)
    print(obj.obj_info)
    from mizore.utilities.Tools import count_operator_nontrivial_terms
    from mizore.method.evolution._utilities import get_first_trotter_gate_count
    print("Nontrivial terms:", count_operator_nontrivial_terms(obj.hamiltonian))
    print("Trotter CNOT:", get_first_trotter_gate_count(obj.hamiltonian))
    #exit(0)
    for quality_cutoff in [0.2]:
        save_path = run_adaptive_product_formula(obj, quality_cutoff=quality_cutoff, n_circuit=10, delta_t=0.2, nM_cutoff=1e11)
        #exit(0)
        benchmark(save_path)
        draw_single(save_path)
