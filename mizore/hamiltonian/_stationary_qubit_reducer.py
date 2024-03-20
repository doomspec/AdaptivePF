from openfermion.ops import QubitOperator
from ..utilities.Tools import get_operator_n_qubit
from ..objective import EnergyObjective
from ..block import HartreeFockInitBlock
from copy import deepcopy

"""
    This file provides the functions to divide or reduce a Hamiltonian
    The main method is 
    get_reduced_energy_obj_with_HF_init(energy_obj:EnergyObjective,location2reduce),
    which can reduce a EnergyObjective by location2reduce. 
    location2reduce is usually set to be the stationary qubits or inactive qubits.
    
    All string_pauli should be like:
    ((0, 'X'), (1, 'X'), (3, 'Y'), (4, 'Y'), (5, 'X'), (6, 'Z'))
"""


def split_pauli_string(string_pauli, location):
    """
    Split a pauli string into two parts
    Return:
    part1_term: correspond to location
    part2_term: other terms
    """
    part1_term = []
    part2_term = list(string_pauli)
    compare_index = 0
    del_num = 0
    for i in range(0, len(location)):
        while compare_index < len(string_pauli):
            if location[i] < string_pauli[compare_index][0]:
                break
            if location[i] > string_pauli[compare_index][0]:
                compare_index += 1
                continue
            if location[i] == string_pauli[compare_index][0]:
                part1_term.append(string_pauli[compare_index])
                part2_term.pop(compare_index - del_num)
                del_num += 1
                break
    return part1_term, part2_term


def get_mapping_by_reduced_location(n_qubit, reduced_location):
    mapping = [-1] * n_qubit
    shift = 0
    for i in range(n_qubit):
        if i in reduced_location:
            shift += 1
        else:
            mapping[i] = i - shift
    return mapping


def get_mapping_by_remain_location(n_qubit, remain_location):
    """
    convert ((1, 'X'), (3, 'Y'),  (5, 'X'), (6, 'Z'))  location = (1,2,3,5,6)to 
    ((0, 'X'), (2, 'Y'),  (3, 'X'), (4, 'Z'))
    """
    if len(remain_location) == 0:
        return []

    mapping = [0] * (remain_location[len(remain_location) - 1] + 1)
    j = 0
    for i in remain_location:
        mapping[i] = j
        j += 1

    return mapping


def map_string(pauli_string, mapping):
    new_string = []
    for term in pauli_string:
        new_string.append((mapping[term[0]], term[1]))
    return new_string


def evaluate_string_list(init_state_string, string_list):
    # inner_prod={"X":[0,0],"Y":[0,0],"Z":[1,-1]}
    res = []
    for string in string_list:
        # print(string)
        coeff = 1
        for term in string:
            if term[1] != "Z":
                coeff = 0
                break
            else:
                if init_state_string[term[0]] == 1:
                    coeff *= -1
        # print(coeff)
        res.append(coeff)
    return res


def get_reduced_energy_obj_with_HF_init(energy_obj: EnergyObjective, location2reduce, location2keep=None,
                                        relabel_qubits=True):
    """
    Args:
        energy_obj: The energy objective to reduce 
        location2reduce: a list like [2,4], which indicates the qubits to reduce 
    """
    if location2keep != None:
        location2reduce = list(set(range(energy_obj.n_qubit)).difference(location2keep))
    location2reduce.sort()
    init_X_qsubset = energy_obj.init_block.qsubset
    init_state_string = [0] * energy_obj.n_qubit
    for i in init_X_qsubset:
        init_state_string[i] = 1
    new_n_qubit = energy_obj.n_qubit - len(location2reduce)
    new_hamiltonian = get_reduced_operator(energy_obj.hamiltonian, location2reduce, init_state_string,
                                           relabel_qubits=relabel_qubits)

    mapping = get_mapping_by_reduced_location(energy_obj.n_qubit, location2reduce)
    new_init_X_qsubset = []
    for i in init_X_qsubset:
        new_index = mapping[i]
        if new_index != -1:
            new_init_X_qsubset.append(new_index)
    return EnergyObjective(new_hamiltonian, new_n_qubit, HartreeFockInitBlock(new_init_X_qsubset),
                           obj_info=energy_obj.obj_info)


def get_reduced_operator(hamiltonian, location2reduce, init_state_string, relabel_qubits=True):
    """
    location2reduce should be sorted
    example: [0,2,3,4]
    """
    n_qubit = get_operator_n_qubit(hamiltonian)
    string_reduce_list = []
    string_remain_list = []
    coff_list = []
    new_hamiltonian = QubitOperator()

    for pauli_and_coff in hamiltonian.get_operators():
        for string_pauli in pauli_and_coff.terms:
            string_reduce, string_remain = split_pauli_string(string_pauli, location2reduce)
            string_reduce_list.append(string_reduce)
            string_remain_list.append(string_remain)
            coff_list.append(pauli_and_coff.terms[string_pauli])

    list_reduced_coeff = evaluate_string_list(init_state_string, string_reduce_list)
    new_string_list = []
    mapping = get_mapping_by_reduced_location(n_qubit, location2reduce)

    if relabel_qubits:
        for i in range(0, len(string_remain_list)):
            new_string_list.append(map_string(string_remain_list[i], mapping))
    else:
        new_string_list = string_remain_list

    for i in range(0, len(coff_list)):
        coff_list[i] = coff_list[i] * list_reduced_coeff[i]
        if string_remain_list[i] == []:
            new_hamiltonian += coff_list[i] * QubitOperator(" ")
        else:
            new_hamiltonian += coff_list[i] * QubitOperator(new_string_list[i])

    return new_hamiltonian
