from .Iterators import iter_qsubset_pauli_of_operator
from openfermion.ops import QubitOperator
import numpy as np
import copy

number2pauli_name = ["I", "X", "Y", "Z"]


def random_list(start, stop, length):
    """
    This function generates a array of random float numbers with certain length 
    The numbers are randomly chosen from start to stop
    """
    import random
    if length >= 0:
        length = int(length)
    start, stop = (start, stop) if start <= stop else (stop, start)
    random_list = []
    for _i in range(length):
        random_list.append(random.uniform(start, stop))
    return random_list


def get_random_coeff_operator(operator, start, end):
    import random
    new_operator = QubitOperator()
    for pauli_and_coeff in operator.get_operators():
        for string_pauli in pauli_and_coeff.terms:
            new_operator += random.uniform(start, end) * QubitOperator(string_pauli)
    return new_operator


def normalize_operator_coeff(operator, average):
    import random
    coeff_sum = 0
    for pauli_and_coeff in operator.get_operators():
        for string_pauli in pauli_and_coeff.terms:
            coeff_sum += abs(pauli_and_coeff.terms[string_pauli])
    new_operator = copy.copy(operator) * (average / coeff_sum)
    return new_operator


def get_operator_qsubset(operator: QubitOperator):
    qubit_set = [False] * 100  # Assume there are at most 100 qubits
    for qsubset, _pauli in iter_qsubset_pauli_of_operator(operator):
        for i in qsubset:
            qubit_set[i] = True
    qsubset = []
    for i in range(len(qubit_set)):
        if qubit_set[i]:
            qsubset.append(i)
    return qsubset


def get_operator_n_qubit(operator: QubitOperator):
    n_qubit = -1
    for pauli_and_coeff in operator.get_operators():
        for string_pauli in pauli_and_coeff.terms:
            if len(string_pauli) == 0:
                continue
            highest_index = string_pauli[len(string_pauli) - 1][0]
            if highest_index > n_qubit:
                n_qubit = highest_index
    return n_qubit + 1


def get_operator_nontrivial_term_weight_sum(operator: QubitOperator):
    coeff_sum = 0
    for pauli_and_coeff in operator.get_operators():
        for string_pauli in pauli_and_coeff.terms:
            if len(string_pauli) == 0:
                print("const term", string_pauli)
                continue
            coeff_sum += abs(pauli_and_coeff.terms[string_pauli])
    return coeff_sum

def count_operator_nontrivial_terms(operator: QubitOperator):
    n_term = 0
    for pauli_and_coeff in operator.get_operators():
        for string_pauli in pauli_and_coeff.terms:
            if len(string_pauli) == 0:
                continue
            n_term+=1
    return n_term


PauliI = np.array([[1, 0], [0, 1]], np.complex)
PauliX = np.array([[0, 1], [1, 0]], np.complex)
PauliY = np.array([[0, -1j], [1j, 0]], np.complex)
PauliZ = np.array([[1, 0], [0, -1]], np.complex)
pauli_dict = {"I": PauliI, "X": PauliX, "Y": PauliY, "Z": PauliZ}


def qubit_operator2matrix(n_qubit, hamiltonian: QubitOperator):
    n_dim = 2 ** n_qubit
    mat = np.zeros((n_dim, n_dim), complex)
    for pauli_and_coff in hamiltonian.get_operators():
        for string_pauli in pauli_and_coff.terms:
            if not string_pauli:
                mat += pauli_and_coff.terms[string_pauli] * np.eye(n_dim)
                continue
            string_pauli_i = 0
            if string_pauli[0][0] == 0:
                pauli0 = pauli_dict[string_pauli[0][1]]
                string_pauli_i += 1
            else:
                pauli0 = PauliI
            for i in range(1, n_qubit):
                if string_pauli_i == len(string_pauli):
                    pauli0 = np.kron(PauliI, pauli0)
                    continue
                if i == string_pauli[string_pauli_i][0]:
                    pauli0 = np.kron(pauli_dict[string_pauli[string_pauli_i][1]], pauli0)
                    string_pauli_i += 1
                else:
                    pauli0 = np.kron(PauliI, pauli0)

            mat += pauli_and_coff.terms[string_pauli] * pauli0
    return mat


def pauliword2string(pauli):
    string = ""
    for i in pauli:
        string += number2pauli_name[i]
    return string
