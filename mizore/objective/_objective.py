class Objective:
    """
    The base class of objectives.
    A objective should include essential information for a problem solver
    Usually, an objective should be able to generate a cost function for the solver to optimize
    """

    def __init__(self, n_qubit=-1, init_block=None):
        self.obj_info = {}
        self.obj_type = ""
        self.init_block = init_block
        self.n_qubit = n_qubit
        return

    def get_cost(self):
        return


class CostFunction:
    def __init__(self):
        return

    def get_cost_obj(self, circuit):
        return

    def get_cost_value(self, circuit):
        return
