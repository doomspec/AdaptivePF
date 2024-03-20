from ._block import Block
from ..utilities.Tools import pauliword2string
from ._compositive_block import CompositiveBlock
from ._block_circuit import BlockCircuit
from ._pauli_gates_block import PauliGatesBlock
from ._minor_blocks import GlobalPhaseBlock
from ._matrix_block import MatrixBlock
from ..backend.gate import PauliRotation
from math import pi


class RotationEntangler(Block):
    """Entangler of the form: e^{iPt}
    Attributes:
        qsubset: The subset of qubits in the wave function that the entangler applies on
        pauliword: The Pauli word P in e^{iPt}
        parameter: t in e^{iPt}
    """
    n_parameter = 1
    IS_INVERSE_DEFINED = True
    IS_DERIVATIVE_DEFINED = True

    derivative_matrix_block=MatrixBlock([0],[[0.5j,0],[0,0.5j]])

    def __init__(self, qsubset, pauliword, init_angle=0.0,temp=1):
        Block.__init__(self, n_parameter=1)
        self.pauliword = pauliword
        self.parameter = [init_angle]
        self.qsubset = qsubset
        self.temp=temp

    def get_forward_gates(self, parameter):
        return [PauliRotation(self.qsubset, self.pauliword, parameter[0] + self.parameter[0])]

    def get_inverse_gates(self, parameter):
        return [PauliRotation(self.qsubset, self.pauliword, -parameter[0] - self.parameter[0])]

    def get_derivative_block(self, para_position):
        bc = BlockCircuit(0)  # Add self as the initial block
        bc.add_block_without_copy(self)
        bc.add_block_without_copy(PauliGatesBlock(self.qsubset, self.pauliword))
        #bc.add_block_without_copy(GlobalPhaseBlock(init_angle=pi / 2))
        #bc.add_block_without_copy(MatrixBlock([0], [[1j, 0], [0, 1j]]))
        bc.add_block_without_copy(MatrixBlock([0],[[1j/2,0],[0,1j/2]])) # TODO
        return CompositiveBlock(bc)

    def get_gate_used(self):
        return {"CNOT": (len(self.qsubset) - 1) * 2,
                "SingleRotation": count_single_gate_for_pauliword(self.pauliword) + 1}

    def get_reconstruct_args(self):
        return [self.qsubset,self.pauliword]

    def __str__(self):
        info = self.basic_info_string()
        info += "; Qsubset:" + str(self.qsubset)
        info += "; Pauli:" + pauliword2string(self.pauliword)
        info += "; Para:" + '%.8f' % self.parameter[0]
        return info


def count_single_gate_for_pauliword(pauliword):
    n_rotation = 0
    for term in pauliword:
        if term != 3:
            n_rotation += 1
    return n_rotation * 2
