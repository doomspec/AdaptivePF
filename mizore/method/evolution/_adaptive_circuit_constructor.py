import numpy.linalg

from mizore.classical.sparse_tools import get_sparse_operator
from ._parameter_evolution_utilities import *
from mizore.circuit_construct.utilities import get_bottom_k_indices
from mizore.block import BlockCircuit
from ._utilities import get_n_measurement_with_var_coeff, get_var_coeff_by_lambda, get_n_measurement
from ...utilities.Tools import get_operator_nontrivial_term_weight_sum


class TimeEvolutionConstructor():

    def __init__(self,n_qubit,block_pool,hamiltonian,quality_cutoff ,task_manager, nM_cutoff):
        self.n_qubit=n_qubit
        self.block_pool = list(block_pool)
        self.quality_cutoff=quality_cutoff
        self.task_manager=task_manager
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

    def construct_circuit(self,circuit:BlockCircuit):

        current_quality=999999
        n_measurement = 1e30
        error_expected = self.quality_cutoff
        # current_A_size=0
        current_pool = [block for block in self.block_pool]
        while (current_quality>self.quality_cutoff) or (n_measurement>self.nM_cutoff):

            parameter_directions = get_parameter_directions(circuit)
            current_state_vec = array(circuit.get_quantum_state().get_vector())
            hamil_state_vec = asarray(dot(self.__np_hamil, current_state_vec))[0]
            vec_C = calc_C_vec(parameter_directions, hamil_state_vec)
            mat_A = calc_A_mat(parameter_directions)
            derivative = calc_derivative(mat_A, vec_C)

            # Break if the objective quality has been achieved
            current_quality = calc_quality(derivative, parameter_directions, hamil_state_vec)

            n_measurement = get_n_measurement(derivative, self.hamiltonian_weight, mat_A, vec_C, current_quality)

            if (current_quality<=self.quality_cutoff) and (n_measurement <= self.nM_cutoff):
                break

            # Send evaluation task to workers
            task_series_id = "Qual" + str(id(self) % 100000)
            task_list = []
            for block in current_pool:
                task_list.append(QualityEvalTask(None, None, None, None, None, block, self.hamiltonian_weight))
            for task in task_list:
                self.task_manager.add_task_to_buffer(
                task, task_series_id=task_series_id)
            self.task_manager.smart_flush(public_resource={"parameter_directions":parameter_directions,"mat_A":mat_A, "vec_C":vec_C, "current_state_vec":current_state_vec, "hamil_state_vec":hamil_state_vec})

            # Receive result
            quality_and_var_coeff_list = self.task_manager.receive_task_result(
                task_series_id=task_series_id, progress_bar=False)
            quality_list = [item[0] for item in quality_and_var_coeff_list]
            nM_list = [get_n_measurement_with_var_coeff(item[1], item[0]) for item in quality_and_var_coeff_list]
            deri_l1_list = [norm(item[2], ord=1) for item in quality_and_var_coeff_list]
            deri_list = [item[2] for item in quality_and_var_coeff_list]
            # Select and add the best block to the circuit
            selected_block = None
            if (current_quality > self.quality_cutoff) or \
                    (n_measurement <= self.nM_cutoff and current_quality <= self.quality_cutoff):
                best_index=get_bottom_k_indices(quality_list,1)[0]
                selected_block = current_pool[best_index]
                circuit.add_active_block(selected_block)
                current_quality = quality_list[best_index]
                n_measurement = nM_list[best_index]
                print("quality_construct", current_quality)
            else:
                best_index = get_bottom_k_indices(deri_l1_list, 1)[0]
                selected_block = current_pool[best_index]
                print("deri_norm_l1", numpy.linalg.norm(deri_list[best_index][1:],1))
                circuit.add_active_block(selected_block)
                current_quality = quality_list[best_index]
                n_measurement = nM_list[best_index]
                print("nM_construct", n_measurement)

                current_quality = quality_list[best_index]
                print("quality_construct", current_quality)

            with open("/home/mh-group/mizore_results/adaptive_evolution/construction.log","a") as f:
                print(current_pool[best_index],quality_list[best_index],file=f)

            current_pool.remove(selected_block)
            #print(len(current_pool))

        return True

from mizore.multiprocess import Task
from mizore.backend import QuantumStateFromVec
from numpy import vstack,hstack

class QualityEvalTask(Task):

    def __init__(self, parameter_directions, mat_A, vec_C, current_state_vec, hamil_state_vec, block, hamiltonian_weight):
        Task.__init__(self)
        self.mat_A=mat_A
        self.vec_C=vec_C
        self.current_state_vec=current_state_vec
        self.hamil_state_vec=hamil_state_vec
        self.parameter_directions=parameter_directions
        self.block=block
        self.hamiltonian_weight = hamiltonian_weight

    def run(self):

        n_new_para=self.block.n_parameter
        n_old_para=len(self.vec_C)
        # state=QuantumStateFromVec(self.current_state_vec)
        # Calculate the new directions
        new_parameter_directions = []
        for i in range(n_new_para):
            state = QuantumStateFromVec(self.current_state_vec)
            derivative_block = self.block.get_derivative_block(i)
            derivative_block | state
            new_parameter_directions.append(state.get_vector()) # TODO
            #print("norm",state.get_squared_norm())
            #derivative_block % state
            #print(state.get_vector()[0:10])

        # Calculate the new mat A
        """
        How new_mat_A is composed
        |------|---|
        |      |   |
        |mat_A |   | <- up_right_mat_A
        |------|---|
        |______|___| <- down_right_mat_A
        """

        up_right_mat_A=zeros((n_old_para,n_new_para),dtype=complex)
        for i in range(n_new_para):
            for j in range(n_old_para):
                up_right_mat_A[j][i]=vdot(self.parameter_directions[j],new_parameter_directions[i])
        up_mat_A=hstack([self.mat_A,up_right_mat_A])

        down_right_mat_A = calc_A_mat(new_parameter_directions)
        down_mat_A=hstack([up_right_mat_A.conjugate().T,down_right_mat_A])

        new_mat_A=vstack([up_mat_A,down_mat_A])

        # Calculate the new vec C
        down_vec_C=zeros((n_new_para,),dtype=complex)
        for i in range(n_new_para):
            down_vec_C[i]=vdot(new_parameter_directions[i], self.hamil_state_vec)
        new_vec_C=hstack([self.vec_C,down_vec_C])
        derivative = calc_derivative(new_mat_A, new_vec_C)
        new_parameter_directions=vstack([self.parameter_directions,new_parameter_directions])
        quality = calc_quality(derivative, new_parameter_directions, self.hamil_state_vec)

        var_coeff = get_var_coeff_by_lambda(derivative, self.hamiltonian_weight, new_mat_A, new_vec_C)

        #return (quality, var_coeff)
        return (quality, var_coeff, derivative)
