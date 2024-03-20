from ._circuit_constructor import CircuitConstructor
from ..block import Block, BlockCircuit
from ..pool import BlockPool
from multiprocessing import Process
from copy import copy, deepcopy
from ..block._utilities import *
from ..objective._objective import Objective
from ..runner import TaskManager, OptimizationTask, GradientTask
from ..optimizer import BasinhoppingOptimizer
from ._result_display import save_construction
import time
import numpy
import math

NOT_DEFINED = 999999


class TrialConstructor(CircuitConstructor):
    """
    Greedy circuit constructor which try all the blocks in self.block_pool in each iteration 
    and add the block that decrease the cost most in the *end* of the self.circuit.

    This strategy is used in the following works
    J. Chem. Theory Comput. 2020, 16, 2
    Nat Commun 10, 3007 (2019)

    Attributes:
        construct_obj: a objective that defines the problem (like EnergyObjective for ground state energy)
        block_pool: The BlockPool used for constructing the circuit, will be gone over at each iteration
        circuit: the circuit kept by the constructor, should be constructed when running
        max_n_iter: Max number of blocks to be added in the circuit
        terminate_cost: the cost where the construction stops
        optimizer: a optimizer for parameter optimization
        task_manager: a TaskManager for parallel run of parameter optimization. 
        If left None, a new task manager uses 4 processes will be created and used 
    """

    gradiant_cutoff = 1e-9

    def __init__(self, construct_obj: Objective, block_pool: BlockPool, max_n_iter=100, gradient_screening_rate=0.05,
                 terminate_cost=-NOT_DEFINED, optimizer=BasinhoppingOptimizer(), global_optimizer=None,
                 no_global_optimization=False, task_manager: TaskManager = None, init_circuit=None,
                 project_name="Untitled", not_save=False):

        CircuitConstructor.__init__(self)
        self.circuit = init_circuit
        if self.circuit == None:
            self.circuit = BlockCircuit(construct_obj.n_qubit)
            self.circuit.add_copied_block(construct_obj.init_block)
        self.max_n_iter = max_n_iter
        self.terminate_cost = terminate_cost
        self.gradient_screening_rate = gradient_screening_rate
        self.block_pool = block_pool
        self.not_save = not_save
        self.n_qubit = construct_obj.n_qubit
        self.cost = construct_obj.get_cost()
        self.id = id(self)
        self.optimizer = optimizer
        self.is_failed = False
        if global_optimizer == None:
            self.global_optimizer = BasinhoppingOptimizer(random_initial=0)
        else:
            self.global_optimizer = global_optimizer
        self.no_global_optimization = no_global_optimization
        self.time_string = time.strftime(
            '%m-%d-%Hh%Mm%Ss', time.localtime(time.time()))
        self.project_name = project_name
        self.save_name = project_name + "_" + self.time_string
        self.trial_circuits = []
        if "terminate_cost" in construct_obj.obj_info.keys():
            self.terminate_cost = construct_obj.obj_info["terminate_cost"]

        self.task_manager = task_manager
        self.task_manager_created = False
        if task_manager == None:
            self.task_manager = TaskManager(n_processor=4)
            self.task_manager_created = True
        return

    CONSTRUCTOR_NAME = "GreedyConstructor"

    def get_current_cost(self):
        return self.cost.get_cost_value(self.circuit)

    def run(self, calc_init=True):

        # Initialization
        print("Here is " + self.CONSTRUCTOR_NAME)
        print("Project Name:", self.project_name)
        print("Block Pool Size:", len(self.block_pool.blocks))
        self.init_cost = self.get_current_cost()
        self.current_cost = self.init_cost
        if self.init_cost < self.terminate_cost:
            print("The objective has already been met! Reture the input circuit.")
            return self.circuit
        self.current_cost = self.init_cost
        self.cost_list.append(self.current_cost)
        print("Initial Cost:", self.init_cost)
        self.start_time_number = time.time()
        self.add_time_point()

        # Iterations
        for i_iter in range(self.max_n_iter):
            print("********The " + str(i_iter + 1) + "th Iteration*********")

            # Try to add a block
            self.update_trial_circuits()
            is_succeed = self.update_one_block()
            is_return = False

            # Post-processing
            if is_succeed:
                # Succeed to add new block
                if self.when_terminate_cost_achieved != -1:
                    if self.task_manager_created:
                        self.task_manager.close()
                    is_return = True
            else:
                # Fail to add new block
                print("Circuit update failed")
                print("Suggestion: 1.Use larger block_pool 2.Ground cost may have achieved")
                print("Final Cost:", self.current_cost)
                print("********Final Circuit********")
                print(self.circuit)
                if self.task_manager_created:
                    self.task_manager.close()
                self.is_failed = True
                is_return = True
            self.add_time_point()
            if not self.not_save:
                save_construction(self, self.save_name)
            if is_return:
                return self.circuit
        print("Circuit Construction ended as it has iterated enough times!")
        print("********Final Circuit********")
        print(self.circuit)
        return self.circuit

    def do_global_optimization(self):

        self.circuit.set_all_block_active()
        task = OptimizationTask(self.circuit, self.global_optimizer, self.cost)
        self.current_cost, parameter = task.run()
        self.circuit.adjust_parameter_on_active_position(parameter)

    def update_one_block(self):
        """Try to add a new block
        Return True is succeed, return False otherwise
        """
        trial_result_list = self.do_trial_on_circuits_by_cost_gradient()
        if len(trial_result_list) != 0:
            self.update_circuit_by_trial_result(trial_result_list)
            print("Block added and shown below, cost now is:",
                  self.current_cost)
            print("********New Circuit********")
            print(self.circuit)
            if not self.no_global_optimization:
                print("Doing global optimization on the new circuit")
                self.do_global_optimization()
                print("Global Optimized Cost:", self.current_cost)

            self.cost_list.append(self.current_cost)

            if not self.not_save:
                save_construction(self, self.save_name)

            print("Distance to target cost:",
                  self.current_cost - self.terminate_cost)
            print("Gate Usage:", self.circuit.get_gate_used())
            print("Cost list:", self.cost_list)
            if self.current_cost <= self.terminate_cost:
                self.when_terminate_cost_achieved = len(
                    self.circuit.block_list)
                print("Target cost achieved by",
                      self.when_terminate_cost_achieved, " blocks!")
                print("Construction process ends!")
            return True
        else:
            print("No trial circuit in the list provides a lower cost")
            return False

    def update_trial_circuits(self, block_pool=None):
        if block_pool == None:
            block_pool = self.block_pool
        self.trial_circuits = []
        for block in block_pool:
            trial_circuit = self.circuit.duplicate()
            trial_circuit.add_copied_block(block)
            trial_circuit.set_only_last_block_active()
            self.trial_circuits.append(trial_circuit)

    def do_trial_on_circuits_by_cost_value(self, trial_circuits=None):
        if trial_circuits == None:
            trial_circuits = self.trial_circuits
        task_series_id = "Single Block Optimize " + str(self.id % 100000)
        trial_result_list = []
        for trial_circuit in trial_circuits:
            task = OptimizationTask(trial_circuit, self.optimizer, None)
            self.task_manager.add_task_to_buffer(task, task_series_id=task_series_id)
        self.task_manager.flush(public_resource={"cost": self.cost})
        res_list = self.task_manager.receive_task_result(
            task_series_id=task_series_id, progress_bar=True)
        for i in range(len(trial_circuits)):
            cost, amp = res_list[i]
            cost_descent = self.current_cost - cost
            # See whether the new entangler decreases the cost
            if cost_descent > self.gradiant_cutoff:
                trial_circuits[i].adjust_parameter_on_active_position(amp)
                trial_result_list.append((cost, trial_circuits[i]))
        return trial_result_list

    def do_trial_on_circuits_by_cost_gradient(self, trial_circuits=None):
        if trial_circuits == None:
            trial_circuits = self.trial_circuits
        if abs(self.gradient_screening_rate - 1) < 1e-5:
            return self.do_trial_on_circuits_by_cost_value()

        task_series_id = "Gradient " + str(self.id % 100000)

        for trial_circuit in trial_circuits:
            task = GradientTask(trial_circuit, None)
            self.task_manager.add_task_to_buffer(task, task_series_id=task_series_id)

        self.task_manager.flush(public_resource={"cost": self.cost})

        res_list = self.task_manager.receive_task_result(task_series_id=task_series_id, progress_bar=True)
        res_list = [numpy.linalg.norm(res) for res in res_list]
        res_list = numpy.array(res_list)
        n_circuit_to_try = math.ceil(
            self.gradient_screening_rate * len(trial_circuits))
        rank_list = (-1 * res_list).argsort()[:n_circuit_to_try]
        good_circuits = []
        for i in rank_list:
            good_circuits.append(trial_circuits[i])
        return self.do_trial_on_circuits_by_cost_value(trial_circuits=good_circuits)

    def update_circuit_by_trial_result(self, trial_result_list):
        lowest_cost = NOT_DEFINED
        lowest_cost_index = -1
        for i in range(len(trial_result_list)):
            if trial_result_list[i][0] < lowest_cost:
                lowest_cost = trial_result_list[i][0]
                lowest_cost_index = i
        best_result = trial_result_list[lowest_cost_index]
        self.current_cost = best_result[0]
        self.circuit = best_result[1]
