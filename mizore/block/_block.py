from copy import deepcopy


class Block():
    """
    Base class of all the blocks
    A block is a piece of parameterized circuit which achieves certain operation.
    The users may override the methods apply_forward_gate() and apply_inverse_gate()
    to define a new block.

    Blocks may be used in a BlockCircuit or a BlockPool.

    Attributes:
        parameter: The list of parameters of the parametric block
        is_inversed: If true, the function apply() will apply apply_inverse_gate() when called
        IS_INVERSE_DEFINED: Set to be true when apply_inverse_gate() is implemented
    """

    IS_INVERSE_DEFINED = False
    IS_DERIVATIVE_DEFINED = False

    def __init__(self, is_inversed=False, n_parameter=-1):
        self.is_inverse = is_inversed
        self.parameter = []
        self.n_parameter = n_parameter

    def get_gates(self, parameter):
        """
        Apply gates parametrized by parameter on the wave function
        """
        if not self.is_inverse:
            return self.get_forward_gates(parameter)
        else:
            return self.get_inverse_gates(parameter)

    def get_forward_gates(self, parameter):
        return []

    def get_inverse_gates(self, parameter):
        # If there is not inverse operation set, use forward operation as inverse.
        assert self.IS_INVERSE_DEFINED
        return

    def adjust_parameter(self, adjuct_list):
        if len(adjuct_list) != len(self.parameter):
            raise Exception("The number of parameter do not match the block!")
        for i in range(len(adjuct_list)):
            self.parameter[i] += adjuct_list[i]

    def basic_info_string(self):
        info = "Type:" + self.__class__.__name__ + \
               "; Para Num:" + str(self.n_parameter)
        if self.is_inverse:
            info += "; INVERSED"
        return info

    def get_gate_used(self):
        return dict()

    def get_derivative_block(self, para_position):
        assert self.IS_DERIVATIVE_DEFINED
        return

    def duplicate(self):
        return deepcopy(self)

    def get_reconstruct_args(self):
        pass

    def __str__(self):
        return self.basic_info_string()

    def __hash__(self):
        return self.__str__().__hash__()

    def __or__(self, state):
        for gate in self.get_forward_gates([0] * self.n_parameter):
            gate.update_quantum_state(state)

    def __mod__(self, state):
        for gate in self.get_inverse_gates([0] * self.n_parameter):
            gate.update_quantum_state(state)

    def __call__(self, para):
        return SmartGateList(self.get_gates(para))

class  SmartGateList():

    def __init__(self,gate_list):
        self.gate_list=gate_list

    def __or__(self, state):
        for gate in self.gate_list:
            gate.update_quantum_state(state)
