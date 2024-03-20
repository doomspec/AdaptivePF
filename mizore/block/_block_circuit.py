import os
import pickle
from copy import deepcopy
from ._block import Block
from ..backend import QuantumCircuit, get_circuit_by_gate_list, QuantumState


class BlockCircuit:
    """
    Circuit consists of blocks. 
    Blocks are listed in self.block_list to form a circuit. 

    And the class provides functions to produce ParametrizedCircuit for parameter optimizers to process. 
    The users can easily fix some parameters while make the others adjustable.

    Especially, in self.active_positions records the indices of blocks in self.block_list
    whose parameter should be adjustable.

    One can make use of duplicate() to duplicate a circuit

    Attributes:
        block_list: The list of blocks contained in the circuit
    """

    def __init__(self, n_qubit, init_block=None):

        self.block_list = []
        self.n_qubit = n_qubit
        self.active_positions = set()
        if init_block is not None:
            self.add_copied_block(init_block)

    """
    Manage Blocks
    """

    def add_copied_block(self, block: Block):
        self.block_list.append(block.duplicate())

    def add_active_block(self, block: Block):
        self.add_copied_block(block)
        self.activate_last_block()

    def activate_last_block(self):
        self.active_positions.add(len(self.block_list) - 1)

    def add_block_without_copy(self, block: Block):
        if block is not None:
            self.block_list.append(block)
            self.active_positions.add(len(self.block_list) - 1)

    def remove_block(self, position):
        self.block_list.pop(position)
        if position in self.active_positions:
            self.active_positions.remove(position)

    def pop_block(self):
        self.remove_block(len(self.block_list) - 1)

    """
    Active position and number of parameters
    """

    def count_n_parameter_by_positions(self, positions):
        n_parameter = 0
        for i in positions:
            n_parameter += self.block_list[i].n_parameter
        return n_parameter

    def count_n_parameter_on_active_positions(self):
        return self.count_n_parameter_by_positions(self.active_positions)

    def count_n_parameter(self):
        positions = list(range(len(self.block_list)))
        return self.count_n_parameter_by_positions(positions)

    def get_parameter_on_positions(self, positions):
        para_list = []
        for position in positions:
            para_list.extend(self.block_list[position].parameter)
        return para_list

    def get_parameter_on_active_position(self):
        return self.get_parameter_on_positions(self.active_positions)

    def set_only_last_block_active(self):
        self.active_positions = set([len(self.block_list) - 1])

    def set_all_block_active(self):
        self.active_positions = set(range(len(self.block_list)))

    def get_active_n_parameter(self):
        return self.count_n_parameter_by_positions(self.active_positions)

    """
    Parameter update
    """

    def adjust_parameter_by_para_position(self, adjust_value, position):
        n_parameter = 0
        block_position = 0
        in_block_position = 0
        for i in range(len(self.block_list)):
            n_parameter += self.block_list[i].n_parameter
            if n_parameter > position:
                block_position = i
                in_block_position = self.block_list[i].n_parameter - \
                                    n_parameter + position
                break
        self.block_list[block_position].parameter[in_block_position] += adjust_value
        return

    def adjust_parameter_by_block_postion(self, adjust_list, position):
        # To be improved
        self.adjust_parameter_by_block_postion_list(adjust_list, (position,))

    def adjust_parameter_by_block_postion_list(self, adjust_list, positions):
        para_index = 0
        for position in positions:
            for in_block_position in range(self.block_list[position].n_parameter):
                self.block_list[position].parameter[in_block_position] += adjust_list[para_index]
                para_index += 1
        return

    def adjust_parameter_on_active_position(self, adjust_list):
        self.adjust_parameter_by_block_postion_list(
            adjust_list, self.active_positions)

    def adjust_all_parameter_by_list(self, adjust_list):

        if len(adjust_list) != self.count_n_parameter():
            raise Exception(
                "The number of parameters provided does not match the circuit!")
        para_index = 0
        for block_position in range(len(self.block_list)):
            n_block_para = self.block_list[block_position].n_parameter
            for in_block_position in range(n_block_para):
                self.block_list[block_position].parameter[in_block_position] += adjust_list[para_index]
                para_index += 1
        return

    """
    Generate ansatz 
    """

    def get_gate_list_by_parameters(self, positions, parameter):
        para_index = 0
        gate_list = []
        for i in range(len(self.block_list)):
            block: Block = self.block_list[i]
            if i in positions:
                gate_list.extend(block.get_gates(parameter[para_index:para_index + block.n_parameter]))
                para_index += block.n_parameter
            else:
                gate_list.extend(block.get_gates([0.0] * block.n_parameter))
        return gate_list

    def get_circuit_by_parameters(self, positions, parameter):
        gate_list = self.get_gate_list_by_parameters(positions, parameter)
        return get_circuit_by_gate_list(self.n_qubit, gate_list)

    def get_ansatz_by_parameters(self, positions):
        def ansatz(parameter):
            state = QuantumState(self.n_qubit)
            circuit: QuantumCircuit = self.get_circuit_by_parameters(positions, parameter)
            circuit.update_quantum_state(state)
            return state
        return ansatz

    def get_gate_list_on_active_parameters(self, parameter):
        return self.get_gate_list_by_parameters(self.active_positions, parameter)

    def get_ansatz_on_active_parameters(self):
        return self.get_ansatz_by_parameters(self.active_positions)

    def get_ansatz_last_block(self):
        positions = [len(self.block_list) - 1]
        return self.get_ansatz_by_parameters(positions)

    def get_ansatz(self):
        positions = list(range(len(self.block_list)))
        return self.get_ansatz_by_parameters(positions)

    def get_fixed_parameter_ansatz(self):
        return self.get_ansatz_by_parameters([])

    def get_quantum_state(self):
        state=QuantumState(self.n_qubit)
        self | state
        return state

    def get_state_vector(self):
        state=QuantumState(self.n_qubit)
        self | state
        return state.get_vector()

    """
    Other functionality
    """

    def get_gate_used(self):
        gate_used = {}
        for block in self.block_list:
            block_gate = block.get_gate_used()
            for key in block_gate.keys():
                if key in gate_used.keys():
                    gate_used[key] += block_gate[key]
                else:
                    gate_used[key] = block_gate[key]
        return gate_used

    def duplicate(self):
        copy_circuit = BlockCircuit(self.n_qubit)
        for block in self.block_list:
            copy_circuit.add_copied_block(block.duplicate())
        copy_circuit.active_positions = self.active_positions
        return copy_circuit

    def __str__(self):
        info = ""
        n_block = len(self.block_list)
        if n_block != 0:
            if n_block > 100:
                return "Huge circuit that contain {} blocks".format(n_block)
            info += "Block Num:" + str(len(self.block_list)) + \
                    "; Qubit Num:" + str(self.n_qubit) + "\n"
            info += "Block list:"
            for block in self.block_list:
                info += "\n" + str(block)
        else:
            info += "\n" + "This is an Empty circuit. Qubit Num:" + \
                    str(self.n_qubit)
        return info

    def __or__(self, other):
        for gate in self.get_gate_list_by_parameters((), ()):
            gate.update_quantum_state(other)
        # TODO compare the efficiency to the following code
        #self.get_circuit_by_parameters((), ()).update_quantum_state(other)

    def save_self_file(self, path, name):
        save_circuit(self, path, name)

    def read_from_file(self, path):
        self = read_circuit(path)


def mkdir(path):
    is_dir_exists = os.path.exists(path)
    if not is_dir_exists:
        os.makedirs(path)
        return True
    else:
        return False


def read_circuit(path):
    with open(path, "rb") as f:
        bc: BlockCircuit = pickle.load(f)
    return bc


def save_circuit(circuit, path, name):
    mkdir(path)
    full_path = path + "/" + name + ".bc"
    with open(full_path, "wb") as f:
        pickle.dump(circuit, f)
