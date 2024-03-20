from ._circuit_constructor import CircuitConstructor
from ..block import BlockCircuit
from ..block_pool import BlockPool
import time
from ._result_display import save_construction
from ..backend import get_circuit_by_gate_list
from ..backend import QuantumState
from multiprocessing import Process, Manager
from .utilities import get_top_k_indices

diff = 0.00001


class SimpleConstructor(CircuitConstructor):

    def __init__(self, n_qubit, cost_wrapper, block_pool: BlockPool, init_state=None, init_block=None,
                 recorder=None,
                 project_name="Untitled"):
        CircuitConstructor.__init__(self)

        if init_state is None:
            self.current_state = QuantumState(n_qubit)
        else:
            self.current_state = init_state

        if init_block is not None:
            init_block | self.current_state

        self.block_pool = list(block_pool.blocks)
        self.n_qubit = n_qubit
        self.cost_wrapper = cost_wrapper

        self.recorder = recorder
        self.circuit = BlockCircuit(n_qubit, init_block=init_block)

        self.max_n_iter = 9999
        self.terminate_cost = -9999

        self.not_save = True
        self.time_string = time.strftime(
            '%m-%d-%Hh%Mm%Ss', time.localtime(time.time()))
        self.project_name = project_name
        self.save_name = project_name + "_" + self.time_string

        return

    def get_current_cost(self):
        return self.cost_wrapper.get_cost_value_by_state(self.current_state)

    def run(self, calc_init=True):
        # Initialization
        print("Here is " + self.CONSTRUCTOR_NAME)
        print("Project Name:", self.project_name)
        print("Block Pool Size:", len(self.block_pool))
        self.init_cost = self.get_current_cost()
        self.current_cost = self.init_cost
        if self.init_cost < self.terminate_cost:
            print("The objective has already been met!")
            return self.circuit
        self.current_cost = self.init_cost
        self.cost_list.append(self.current_cost)
        print("Initial Cost:", self.init_cost)
        self.start_time_number = time.time()
        self.add_time_point()

        # Iterations
        for i_iter in range(self.max_n_iter):
            print("*" * 10 + "The " + str(i_iter + 1) + "th Iteration" + "*" * 10)

            is_succeed = self.update_state()
            is_return = False

            # Post-processing
            if is_succeed:
                # Succeed to add new block
                if self.when_terminate_cost_achieved != -1:
                    is_return = True
            else:
                # Fail to add new block
                print("Circuit update failed")
                print("Suggestion: 1.Use larger block_pool 2.Ground cost may have achieved")
                print("Final Cost:", self.current_cost)
                print("*" * 10 + "Final Circuit" + "*" * 10)
                print(self.circuit)
                is_return = True

            self.add_time_point()

            if not self.not_save:
                save_construction(self, self.save_name)
            if is_return:
                return self.circuit

        print("Circuit Construction ended as it has iterated enough times!")
        print("*" * 10 + "Final Circuit" + "*" * 10)
        print(self.circuit)
        return self.circuit

    def update_state(self):
        state = self.current_state.copy()


        state_vec=state.get_vector()
        hamil_str=str(self.cost_wrapper.hamiltonian)
        tasks=[block.get_reconstruct_args() for block in self.block_pool]
        from ..multiprocess.utilities import get_task_packages
        task_packages=get_task_packages(tasks)

        gradient_lists = []#[get_new_block_gradients0((state_vec,package,hamil_str)) for package in task_packages]


        from multiprocessing import Pool
        pool = Pool(processes=1)
        print("Hello")
        gradient_lists0 = pool.map(get_new_block_gradients0, [(state_vec,package,hamil_str) for package in task_packages])
        print(gradient_lists0)

        gradients = []
        for gradient_list in gradient_lists:
            gradients.extend(gradient_list)
        print(gradients)

        return False


        scores = [abs(gradient) for gradient in gradients]
        top_k = get_top_k_indices(scores, 1)
        new_block = self.block_pool[top_k[0]].duplicate()
        new_block.parameter = [- gradients[top_k[0]] * 0.1]
        new_block | self.current_state
        self.circuit.add_block_without_copy(new_block)
        return True


def get_new_block_gradients(init_state: QuantumState, block_pool, cost_wrapper):
    init_cost = cost_wrapper.get_cost_value_by_state(init_state)
    print(init_cost)
    diff = 0.00001
    gradients = [0] * len(block_pool)
    i = 0
    for block in block_pool:
        block.parameter = [diff]
        block | init_state
        gradient = (cost_wrapper.get_cost_value_by_state(init_state) - init_cost) / diff
        block % init_state
        block.parameter = [0]
        gradients[i] = gradient
        i += 1
    return gradients

from math import log2
from ..block import RotationEntangler
from qulacs.quantum_operator import create_quantum_operator_from_openfermion_text

def get_new_block_gradients0(input_tuple):
    init_state_vec=input_tuple[0]
    block_pool_args=input_tuple[1]
    openfermion_str=input_tuple[2]
    n_qubit=int(log2(len(init_state_vec)))
    print("Hello1")
    init_state = QuantumState(n_qubit)
    print("Hello2")
    init_state.load(init_state_vec)
    print("Hello3")
    qulacs_hamiltonian = create_quantum_operator_from_openfermion_text(openfermion_str)
    print("Hello4")
    init_cost = qulacs_hamiltonian.get_expectation_value(init_state)
    print("Hello5")
    diff = 0.00001
    gradients = [0] * len(block_pool_args)
    i = 0
    for arg in block_pool_args:
        block=RotationEntangler(*arg)
        block.parameter = [diff]
        block | init_state
        gradient = (qulacs_hamiltonian.get_expectation_value(init_state) - init_cost) / diff
        block % init_state
        block.parameter = [0]
        gradients[i] = gradient
        i += 1
    return gradients
