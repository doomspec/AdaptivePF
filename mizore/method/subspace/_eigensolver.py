
from openfermion import QubitOperator
from mizore.classical import get_operator_matrix
import numpy as np
from numpy import vdot,conjugate,dot,asarray
from scipy.linalg import eigh

class SubspaceEigensolver:

    def __init__(self,vector_list=None,hamiltonian_mat=None):
        self.vector_list=vector_list
        self.hamiltonian_mat=hamiltonian_mat
        self.hamiltonian=None
        self.S_mat = None
        self.H_mat = None
        self.eigvals = None
        self.eigvecs = None
        self.ground_energy = None
        self.ground_state = None
        pass


    def load_block_circuits(self,block_circuits):
        self.vector_list=[bc.get_state_vector() for bc in block_circuits]
        pass

    def load_quantum_states(self,quantum_states):
        self.vector_list=[state.get_vector() for state in quantum_states]
        pass

    def load_hamiltonian_ops(self,hamiltonian_ops:QubitOperator):
        self.hamiltonian_mat=get_operator_matrix(hamiltonian_ops)
        pass


    def calc_S_mat(self):
        n_basis=len(self.vector_list)
        S_mat = np.array([[0.0] * n_basis] * n_basis, dtype=complex)
        for i in range(n_basis):
            for j in range(i, n_basis):
                if i == j:
                    S_mat[i][j] = 1.0
                    continue
                s = vdot(self.vector_list[i],self.vector_list[j])
                S_mat[i][j] = s
                S_mat[j][i] = conjugate(s)
        self.S_mat=S_mat


    def calc_H_mat(self):
        n_basis = len(self.vector_list)
        hamil_dot_vec_list=[None]*n_basis
        for i in range(n_basis):
            hamil_dot_vec_list[i]=asarray(dot(self.hamiltonian_mat,self.vector_list[i]))[0]
        H_mat = np.array([[0.0] * n_basis] * n_basis, dtype=complex)
        for i in range(n_basis):
            for j in range(i, n_basis):
                s = vdot(self.vector_list[i], hamil_dot_vec_list[j])
                H_mat[i][j] = s
                H_mat[j][i] = conjugate(s)
        self.H_mat=H_mat


    def solve(self):
        self.calc_S_mat()
        revise_little_negative(self.S_mat)
        self.calc_H_mat()
        # print(self.H_mat)
        # print(self.S_mat)
        self.eigvals, self.eigvecs = eigh(
            self.H_mat, self.S_mat, eigvals_only=False)
        self.ground_energy = self.eigvals[0]
        self.ground_state = self.eigvecs[:, 0]

        return self.ground_energy



def revise_little_negative(S_mat: np.array):
    eps = 1e-9
    eigv = np.linalg.eigvalsh(S_mat)
    assert eigv[0] > -eps
    if eigv[0] < eps:
        S_mat += np.eye(len(S_mat)) * eps