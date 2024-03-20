from ._block import Block
from qulacs.gate import Pauli


class PauliGatesBlock(Block):
    n_parameter = 0
    IS_INVERSE_DEFINED = True

    def __init__(self, qsubset, pauliword):
        self.pauliword = pauliword
        self.qsubset = qsubset

        Block.__init__(self, n_parameter=0)

    def get_forward_gates(self, parameter):
        return [Pauli(self.qsubset, self.pauliword)]

    def get_inverse_gates(self, parameter):
        return [Pauli(self.qsubset, self.pauliword)]

    def get_gate_used(self):
        return {"SingleRotation": len(self.paulistring)}

    def get_reconstruct_args(self):
        return [self.qsubset,self.pauliword]

    def __str__(self):
        info = self.basic_info_string()
        info += "; PauliString:" + str(self.pauliword)
        return info
