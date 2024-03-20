from ._block import Block
from mizore.backend.gate._state_preparation_gate import StatePreparationGate


class StatePreparationBlock(Block):

    n_parameter = 0
    IS_INVERSE_DEFINED = False

    def __init__(self, state_vector):
        Block.__init__(self, n_parameter=0)
        self.state_preparation_gate = StatePreparationGate(state_vector)

    def get_forward_gates(self, parameter):
        return [self.state_preparation_gate]

    def get_inverse_gates(self, parameter):
        assert False

    def get_gate_used(self):
        return dict()

    def __str__(self):
        info = self.basic_info_string()
        info += "; State reset"
        return info
