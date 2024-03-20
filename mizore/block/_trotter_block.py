from ._block import Block
from openfermion.ops import QubitOperator
from ..utilities.Iterators import iter_coeff_qsubset_pauli_of_operator
from ._rotation_entangler import count_single_gate_for_pauliword
from ..backend.gate import PauliRotation

class TrotterEvolutionBlock(Block):
    n_parameter = 1
    IS_INVERSE_DEFINED = True

    def __init__(self, hamiltonian: QubitOperator, n_trotter_step=1, evolution_time=0):
        Block.__init__(self, n_parameter=1)
        self.n_trotter_step = n_trotter_step
        self.qsubset_pauliword_list = []
        for term in iter_coeff_qsubset_pauli_of_operator(hamiltonian):
            self.qsubset_pauliword_list.append(term)
        self.parameter = [evolution_time]

    def iterate_operator(self, total_evolution_time, is_inverse=False):
        step_evolution_angle = -2*(total_evolution_time) / self.n_trotter_step
        term_iterator = range(len(self.qsubset_pauliword_list))

        if is_inverse:
            term_iterator = list(reversed(term_iterator))
            step_evolution_angle = -step_evolution_angle

        for _step in range(self.n_trotter_step):
            for i in term_iterator:
                coeff, qsubset, pauliword = self.qsubset_pauliword_list[i]
                yield qsubset, pauliword, coeff * step_evolution_angle

    def get_forward_gates(self, parameter):
        gate_list=[PauliRotation(*pauli_rotation) for pauli_rotation in
                   self.iterate_operator(parameter[0] + self.parameter[0])]
        return gate_list

    def get_inverse_gates(self, parameter):
        gate_list = [PauliRotation(*pauli_rotation) for pauli_rotation in
                     self.iterate_operator(parameter[0] + self.parameter[0], is_inverse=True)]
        return gate_list

    def get_gate_used(self):
        n_rotation = 0
        n_CNOT = 0
        for _ecoeff, qsubset, pauliword in self.qsubset_pauliword_list:
            n_rotation += count_single_gate_for_pauliword(pauliword)
            n_CNOT += 2 * (len(qsubset) - 1)
        n_rotation += 1
        return {"CNOT": n_CNOT * self.n_trotter_step, "SingleRotation": n_rotation * self.n_trotter_step}

    def __str__(self):
        info = self.basic_info_string()
        info += "; Trotter_TimeEvolution: T=" + \
                str(self.parameter[0]) + " N_step:" + str(self.n_trotter_step)
        return info
