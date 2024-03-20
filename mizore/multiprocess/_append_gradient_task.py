from ._task import Task
from ..block import BlockCircuit
from ..backend import get_backend_operator, QuantumState

class AppendGradientTask(Task):

    def __init__(self, circuit: BlockCircuit, blocks2append, hamiltonian, diff=1e-6):
        Task.__init__(self)
        self.circuit = circuit
        self.blocks2append = blocks2append
        self.hamiltonian = hamiltonian
        self.diff = diff

    def run(self):
        old_energy = 0
        gradient_list = []
        backend_hamil = get_backend_operator(self.hamiltonian)
        # Initialize the wavefunction
        state = QuantumState(self.circuit.n_qubit)
        self.circuit | state

        old_energy = backend_hamil.get_expectation_value(state)

        for block in self.blocks2append:
            gradient = []
            for i in range(block.n_parameter):
                diff_para = [0.0] * block.n_parameter
                diff_para[i] = self.diff
                block(diff_para) | state
                energy = backend_hamil.get_expectation_value(state).real
                gradient.append((energy - old_energy) / self.diff)
                block.is_inverse=True
                block(diff_para) | state
                block.is_inverse=False

            gradient_list.append(gradient)

        return gradient_list
