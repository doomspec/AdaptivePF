from sys import path
path.append("/home/mh-group/zijian/Remote/Mizore/")

from mizore.method.evolution import run_adaptive_product_formula
from benchmark import benchmark
from mizore.hamiltonian import make_example_H4
from draw_single import draw_single



if __name__ == '__main__':
    root_path = "/home/mh-group/mizore_results/adaptive_evolution/h4_new"

    for quality_cutoff in [0.2]:
        obj = make_example_H4()
        coeff = [value for value in obj.hamiltonian.terms.values()]

        from mizore.utilities.Tools import count_operator_nontrivial_terms
        from mizore.method.evolution._utilities import get_first_trotter_gate_count

        print("Nontrivial terms:", count_operator_nontrivial_terms(obj.hamiltonian))
        print("Trotter CNOT:", get_first_trotter_gate_count(obj.hamiltonian))

        save_path = run_adaptive_product_formula(obj, quality_cutoff=quality_cutoff, n_circuit=10, delta_t=0.2, nM_cutoff=1e11)
        #exit(0)
        benchmark(save_path, trotter=[])
        draw_single(save_path)
