import time
import json
import pickle
import os
from ._parameter_evolver import TimeEvolutionEvolver
from ._adaptive_circuit_constructor import TimeEvolutionConstructor
import numpy as np
from mizore.block import BlockCircuit


class QualityEnsuredEvolver():

    default_para_dict = {"project_name": "Untitled", "is_to_save": False,
                         "n_circuit": 10, "delta_t": 0.1, "special_save_name": None,"run_info_more":None, "save_root":"~/"}

    def __init__(self, init_circuit: BlockCircuit, constructor:TimeEvolutionConstructor,
                 evolver: TimeEvolutionEvolver, task_manager=None, para_dict=None):

        if para_dict is None:
            para_dict = dict()

        _para_dict = QualityEnsuredEvolver.default_para_dict.copy()
        _para_dict.update(para_dict)

        self.para_dict = _para_dict

        self.n_qubit = init_circuit.n_qubit
        self.circuit = init_circuit
        self.evolver = evolver
        self.constructor=constructor

        self.task_manager = task_manager

        self.project_name = _para_dict["project_name"]
        self.delta_t = _para_dict["delta_t"]
        self.n_circuit = _para_dict["n_circuit"]
        self.is_to_save = _para_dict["is_to_save"]
        self.special_save_name = _para_dict["special_save_name"]
        self.run_info_more = _para_dict["run_info_more"]
        self.save_root=_para_dict["save_root"]
        self.time_string = time.strftime(
            '%m-%d-%Hh%Mm%Ss', time.localtime(time.time()))

        self.quality_list = []
        self.evolution_time_list = []
        self.n_block_change = []
        self.total_time_evolved = 0

        self.circuit_list = [self.circuit.duplicate()]
        self.total_time_evolved = 0

        if self.special_save_name is None:
            self.save_name = self.project_name+"_"+self.time_string
        else:
            self.save_name = self.special_save_name

        self.save_path = os.path.join(_para_dict["save_root"],self.save_name)
        self.save_self_info()

        if self.is_to_save:
            self.circuit.save_self_file(self.save_path+"/circuits", "0")

    """
    Main module
    """

    def run(self, construct_first=False, progress_bar=True):

        last_evolved_time = -1
        total_time_evolved_list = []
        construct_needed = construct_first  # Default to be True

        while True:

            if construct_needed:
                is_success = self.constructor.construct_circuit(self.circuit)
                if is_success is not True:
                    print("Circuit constructor can not reach the object quality, consider using a larger pool!")
                    assert False
                self.n_block_change.append((self.total_time_evolved, len(
                    self.circuit.block_list), self.circuit.get_gate_used()))

            local_time_to_evolve = self.delta_t - \
                (self.total_time_evolved % self.delta_t)
            if local_time_to_evolve < 1e-12:
                local_time_to_evolve = self.delta_t

            print("local_time_to_evolve", local_time_to_evolve)

            # Evolve the circuit
            # The Evolver will change the circuit directly
            evolved_time = self.evolver.evolve_circuit(self.circuit,local_time_to_evolve)

            self.quality_list = np.append(
                self.quality_list, self.evolver.quality_list)
            self.evolution_time_list = np.append(
                self.evolution_time_list, self.total_time_evolved+self.evolver.evolution_time_list)

            self.total_time_evolved += evolved_time

            print("Time evolved:", evolved_time)

            construct_needed = (
                abs(evolved_time-local_time_to_evolve) >= 1e-10)

            # Save the circuit
            if (not construct_needed) and (self.total_time_evolved > last_evolved_time):

                self.circuit_list.append(self.circuit.duplicate())
                last_evolved_time = self.total_time_evolved
                total_time_evolved_list.append(self.total_time_evolved)
                print("Circuit added, progress:{}/{}".format(len(self.circuit_list)-1,self.n_circuit))

                if self.is_to_save:
                    self.circuit.save_self_file(
                        self.save_path+"/circuits", str(len(self.circuit_list)-1))
                    self.save_run_status_info()

            # Exit when enough circuits are generated
            if len(self.circuit_list) == self.n_circuit+1:
                return self.circuit_list

    """
    Save and read of a run
    """

    def save_self_info(self):
        """
        Generate a log as a JSON file
        """
        mkdir(self.save_path)

        para_dict = {}
        for key in QualityEnsuredEvolver.default_para_dict.keys():
            para_dict[key] = self.__dict__[key]
        path = self.save_path + "/para_dict.json"
        with open(path, "w") as f:
            json.dump(para_dict, f)
        path = self.save_path + "/hamiltonian.pickle"
        with open(path, "wb") as f:
            pickle.dump(self.evolver.hamiltonian, f)
        path = self.save_path + "/pool.pickle"
        with open(path, "wb") as f:
            pickle.dump(self.constructor.block_pool, f)

    def get_run_status_info_dict(self):
        log_dict = {"quality_cutoff":self.evolver.quality_cutoff,"hamiltonian_weight":self.evolver.hamiltonian_weight,"n_block_change": self.n_block_change,"quality_list": self.quality_list.tolist(),
                    "evolution_time_list": self.evolution_time_list.tolist()}
        return log_dict

    def save_run_status_info(self):
        path = self.save_path + "/run_status_info.json"
        with open(path, "w") as f:
            json.dump(self.get_run_status_info_dict(), f)


def mkdir(path):
    is_dir_exists = os.path.exists(path)
    if not is_dir_exists:
        os.makedirs(path)
        return True
    else:
        return False
