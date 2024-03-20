from ._block import Block
from ..utilities.Operations import apply_time_evolution
from ..utilities.Tools import get_operator_qsubset
from ..classical import get_operator_matrix
from numpy.linalg import eigh
from numpy import dot,diag,exp,array
from copy import copy
from ..backend.gate import DenseMatrix

class TimeEvolutionBlock(Block):
    n_parameter = 1
    IS_INVERSE_DEFINED = True

    def __init__(self, hamiltonian, init_angle=0):

        self.qsubset = get_operator_qsubset(hamiltonian)
        Block.__init__(self, n_parameter=1)
        self.n_qubit=len(self.qsubset)
        self.hamiltonian = hamiltonian
        self.parameter = [init_angle]
        # We will diagonalize the Hamiltonian H into H=P*D*P^dagger
        self.vec_D=None
        self.mat_P=None

    def diagonalize(self):
        hamiltonian_mat=get_operator_matrix(self.hamiltonian,n_qubits=self.n_qubit)
        self.vec_D, self.mat_P = eigh(hamiltonian_mat)

    def get_evolution_operator(self,evolve_time):
        # This part may change to use Sparse Matrix
        time_evol_op = dot(dot(self.mat_P, diag(exp(-1j * evolve_time * self.vec_D))),
                              self.mat_P.T.conj())
        return time_evol_op

    def get_forward_gates(self, parameter):
        if self.vec_D is None:
            self.diagonalize()
        time_evol_op=self.get_evolution_operator(parameter[0]+self.parameter[0])
        time_evol_gate = DenseMatrix(self.qsubset, time_evol_op)
        return [time_evol_gate]

    def get_inverse_gates(self, parameter):
        if self.vec_D is None:
            self.diagonalize()
        time_evol_op = self.get_evolution_operator(-(parameter[0]+self.parameter[0]))
        time_evol_gate = DenseMatrix(self.qsubset, time_evol_op)
        return [time_evol_gate]

    def duplicate(self):
        new_block=TimeEvolutionBlock(self.hamiltonian,init_angle=self.parameter[0])
        new_block.is_inverse=self.is_inverse
        return new_block

    def get_gate_used(self):
        return {"TimeEvolution": 1}

    def __str__(self):
        info = self.basic_info_string()
        info += "; TimeEvolution: T=" + str(self.parameter[0])
        return info
