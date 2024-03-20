from . import BlockPool
from ..block import RotationEntangler
from ..utilities.Iterators import iter_qsubset_odd_Y_pauli_by_length, iter_all_qsubset_pauli_by_length, \
    iter_qsubset_pauli_of_operator
from copy import copy


def all_rotation_pool(n_qubit, max_length=-1, only_odd_Y_operators=True):
    """
    Block block_pool contains all the possible RotationEntangler.
    Args:
        only_odd_Y_operators: If true, only RotationEntangler whose Pauli word only has odd number of Y operators will be added.
        Even Y entanglers are excluded because they commute with Hamiltonian without imaginary terms (usually because of the absense of magnetic field), which means
        d/dt <psi|e^{-iPt}He^{iPt}|psi> = 0,
        for all the state |psi>
    """
    if max_length == -1:
        max_length = n_qubit

    if only_odd_Y_operators:
        for length in [i + 1 for i in range(max_length)]:
            for qsubset, pauli in iter_qsubset_odd_Y_pauli_by_length(length, range(n_qubit)):
                yield RotationEntangler(qsubset, pauli)
    else:
        for length in [i + 1 for i in range(max_length)]:
            for qsubset, pauli in iter_all_qsubset_pauli_by_length(length, range(n_qubit)):
                yield RotationEntangler(qsubset, pauli)


def make_pauli_imaginary(pauli):
    """
    Modify the Pauli word by replacing one Y by X or one X by Y.
    Return the modified Pauli word.
    """
    new_pauli = copy(pauli)
    for i in range(len(pauli)):
        if new_pauli[i] == 1:
            new_pauli[i] = 2
            return new_pauli
        if new_pauli[i] == 2:
            new_pauli[i] = 1
            return new_pauli
    return None


def quasi_imaginary_evolution_rotation_pool(hamiltonian):
    """
    Entangler block_pool inspired by graph from arXiv:1908.09533v1.
    The block_pool consists of the RotationEntangler with Pauli words modified from the Hamiltonian.
    The modification of Pauli words is to replace one Y by X or one X by Y.
    """
    # print(hamiltonian)

    for qsubset, pauli in iter_qsubset_pauli_of_operator(hamiltonian):
        pauli = make_pauli_imaginary(pauli)
        if pauli != None:
            yield RotationEntangler(qsubset, pauli)


def hamiltonian_rotations_pool(hamiltonian):
    for qsubset, pauli in iter_qsubset_pauli_of_operator(hamiltonian):
        yield RotationEntangler(qsubset, pauli)
