from . import Block
from ..backend.gate import GlobalPhaseGate
from math import pi


class GlobalPhaseBlock(Block):
    n_parameter = 1
    IS_INVERSE_DEFINED = True
    IS_LOCALIZE_AVAILABLE = True
    IS_DERIVATIVE_DEFINE = True

    def __init__(self, init_angle=0.0):
        Block.__init__(self, n_parameter=1)
        self.parameter = [init_angle]

    def get_forward_gates(self, parameter):
        return [GlobalPhaseGate(parameter[0] + self.parameter[0])]

    def get_inverse_gates(self, parameter):
        return [GlobalPhaseGate(-parameter[0] - self.parameter[0])]

    def get_derivative_block(self, para_position):
        return GlobalPhaseBlock(init_angle=self.parameter[0] + pi / 2)

    def get_gate_used(self):
        return {}

    def __str__(self):
        info = self.basic_info_string()
        info += "; T=" + str(self.parameter[0])
        return info