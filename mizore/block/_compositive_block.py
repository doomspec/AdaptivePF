from ._block import Block
from ._block_circuit import BlockCircuit
from ._utilities import get_inverse_circuit


class CompositiveBlock(Block):
    """
    To define a Block by a block circuit
    """
    IS_INVERSE_DEFINED = True

    def __init__(self, circuit: BlockCircuit):
        Block.__init__(self, n_parameter=0)
        self.circuit = circuit
        self.inverse_circuit = get_inverse_circuit(circuit)
        self.n_parameter = self.circuit.count_n_parameter()
        self.parameter = [0.0] * self.n_parameter

    def get_forward_gates(self, parameter):
        return self.circuit.get_gate_list_on_active_parameters(parameter)

    def get_inverse_gates(self, parameter):
        return self.inverse_circuit.get_gate_list_on_active_parameters([-para for para in parameter])

    def __str__(self):
        info = self.basic_info_string() + "\n"
        info += "start \n"
        for block in self.circuit.block_list:
            info += str(block) + "\n"
        info += "end"
        return info
