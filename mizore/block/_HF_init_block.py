from ._block import Block
from qulacs.gate import Pauli
import json


class HartreeFockInitBlock(Block):
    """Apply X gates on a few of qubits. Usually for getting Hartree-Fock qubit wavefunction
    Attributes:
        qsubset: should be the qubits where X gates to be applied
    """
    n_parameter = 0
    IS_INVERSE_DEFINED = True

    def __init__(self, qsubset):
        Block.__init__(self, n_parameter=0)
        self.qsubset = qsubset

    def get_forward_gates(self, parameter):
        return [Pauli(self.qsubset, [1] * len(self.qsubset))]

    def get_inverse_gates(self, parameter):
        return [Pauli(self.qsubset, [1] * len(self.qsubset))]

    def get_gate_used(self):
        return {"SingleRotation": len(self.qsubset)}

    def get_reconstruct_args(self):
        return [self.qsubset]

    def __str__(self):
        info = self.basic_info_string()
        info += "; Qsubset:" + str(self.qsubset)
        return info
