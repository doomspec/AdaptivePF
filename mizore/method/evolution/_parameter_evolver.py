from mizore.multiprocess import TaskManager
from mizore.block import BlockCircuit
from numpy.linalg import norm
from numpy import array
from ._parameter_evolution_utilities import calc_derivative_quality
from mizore.classical.sparse_tools import get_sparse_operator
from mizore.hamiltonian._utilities import get_operator_nontrivial_term_weight_sum
import numpy as np
import math

from ._utilities import get_n_measurement


class TimeEvolutionEvolver():

    def __init__(self, n_qubit, hamiltonian, quality_cutoff=0.001, diff=1e-4, stepsize=1e-3,
                 inverse_evolution=False, task_manager: TaskManager = None, verbose=False,
                 fig_path=None, nM_cutoff=5e8):

        self.n_qubit = n_qubit
        self.stepsize = stepsize
        self.diff = diff
        self.inverse_evolution = inverse_evolution
        self.verbose = verbose
        self.fig_path = fig_path
        self.task_manager = task_manager

        self.quality_cutoff = quality_cutoff
        self.evolution_time_list = None
        self.quality_list = None

        self.__hamil = hamiltonian
        self.__np_hamil = get_sparse_operator(hamiltonian, n_qubits=self.n_qubit).todense()

        self.hamiltonian_weight = get_operator_nontrivial_term_weight_sum(hamiltonian)
        self.nM_cutoff = nM_cutoff

    @property
    def hamiltonian(self):
        return self.__hamil

    @hamiltonian.setter
    def hamiltonian(self, hamiltonian):
        self.__hamil = hamiltonian
        self.__np_hamil = get_sparse_operator(hamiltonian, n_qubits=self.n_qubit).todense()

    def evolve_circuit(self, circuit: BlockCircuit, evolution_time, max_n_step=1000):

        n_parameter = circuit.count_n_parameter_on_active_positions()
        self.evolution_time_list = []
        self.quality_list = []

        if n_parameter == 0:
            print("ATTENTION: No active block in the circuit")
            return 0

        evolved_time = 0
        n_step = 0

        while True:
            derivative, quality, mat_A, vec_C = calc_derivative_quality(self.__np_hamil, circuit)
            error_expected = self.quality_cutoff
            nM = get_n_measurement(derivative, self.hamiltonian_weight, mat_A, vec_C, quality)

            print("nM_evolving", nM, "cutoff", self.nM_cutoff)
            print("deri_norm_ratio", norm(derivative[1:],ord=1)/self.hamiltonian_weight)

            if nM > self.nM_cutoff:
                break

            # Break when quality is above cutoff
            if quality >= self.quality_cutoff:
                break

            self.evolution_time_list.append(evolved_time)
            self.quality_list.append(quality)

            delta_t_evolve = self.stepsize

            if evolved_time + delta_t_evolve >= evolution_time:
                delta_t_evolve = evolution_time - evolved_time

            para_shift = derivative * delta_t_evolve
            evolved_time += delta_t_evolve

            n_step += 1
            if self.inverse_evolution:
                para_shift = -1 * para_shift

            circuit.adjust_parameter_on_active_position(para_shift)

            if evolved_time >= evolution_time - 1e-10:
                break
            if n_step >= max_n_step:
                print("ATTENTION: n_step>=max_n_step ! time evolved in this step:", evolved_time)

        self.evolution_time_list = array(self.evolution_time_list)
        self.quality_list = array(self.quality_list)

        return evolved_time


