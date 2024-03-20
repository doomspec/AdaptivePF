from multiprocessing import Process
import time

NOT_DEFINED = 999999


class CircuitConstructor():
    """
    The base class of all circuit constructor
    See GreedyConstructor to learn how to use
    """

    CONSTRUCTOR_NAME = "CircuitConstructor"

    def __init__(self):
        self.block_pool = None
        self.circuit = None  # Should be a BlockCircuit
        self.init_cost = NOT_DEFINED
        self.terminate_cost = -NOT_DEFINED
        self.max_n_block = NOT_DEFINED
        self.when_terminate_cost_achieved = -1
        self.current_cost = NOT_DEFINED
        self.init_operator = None
        self.start_time_number = NOT_DEFINED
        self.important_log_list = []
        self.run_time_list = []
        self.cost_list = []
        return

    def add_time_point(self):
        self.run_time_list.append(time.time() - self.start_time_number)

    def run(self):
        pass
