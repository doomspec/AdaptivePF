from sys import path
path.append("/home/mh-group/zijian/Remote/Mizore/")

from mizore.method.evolution import run_adaptive_product_formula
from mizore.hamiltonian.random_coeff import get_random_coeff_operator, normalize_operator_coeff
from mizore.hamiltonian.lattice.ising import full_connected_transverse_field_ising

from benchmark import benchmark
from draw_single import draw_single

if __name__ == '__main__':

    for i in range(1):

        n_qubit = 12
        obj = full_connected_transverse_field_ising(n_qubit, 1, 1)
        random_hamil = get_random_coeff_operator(obj.hamiltonian, -1, 1)
        normalized_hamil = normalize_operator_coeff(random_hamil,
                                                    (((n_qubit - 1) * n_qubit) // 2 + n_qubit) // 2)  # hamil_weight=n_qubit
        obj.hamiltonian = normalized_hamil

        from mizore.utilities.Tools import count_operator_nontrivial_terms
        from mizore.method.evolution._utilities import get_first_trotter_gate_count

        print("Nontrivial terms:", count_operator_nontrivial_terms(obj.hamiltonian))
        print("Trotter CNOT:", get_first_trotter_gate_count(obj.hamiltonian))



        save_path = run_adaptive_product_formula(obj, quality_cutoff=0.2, n_circuit=10, delta_t=0.1, n_worker=70)
        exit(0)
        benchmark(save_path)
        draw_single(save_path)
