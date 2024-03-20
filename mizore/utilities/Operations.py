from projectq.ops import CZ, H, X, Y, Z, All, Measure, CNOT, Z, Rz, Ry, Rx, C, TimeEvolution
import math
import projectq
from openfermion.ops import QubitOperator

PAULI_CHAR2OPERATION = {"X": X, "Y": Y, "Z": Z}

"""
The common quantum operations that are used in Mizore
"""


def apply_time_evolution(hamiltonian: QubitOperator, time, wavefunction):
    projectq_qubit_operator = projectq.ops.QubitOperator()
    for term, coefficient in hamiltonian.terms.items():
        projectq_qubit_operator.terms[term] = coefficient
    # print(time)
    TimeEvolution(time, projectq_qubit_operator) | wavefunction


def apply_X_gates(qsubset, wavefunction):
    for i in qsubset:
        X | wavefunction[i]


def apply_H_gates(qsubset, wavefunction):
    for i in qsubset:
        H | wavefunction[i]


def apply_CZ_gates(pairset, wavefunction):
    for pair in pairset:
        CZ | (wavefunction[pair[0]], wavefunction[pair[1]])


def apply_Pauli_gates(paulistring, wavefunction):
    for term in paulistring:
        PAULI_CHAR2OPERATION[term[1]] | wavefunction[term[0]]


def CNOT_entangler(wavefunction, qsubset):
    for i in range(len(qsubset) - 1):
        CNOT | (wavefunction[qsubset[i]], wavefunction[qsubset[i + 1]])
    CNOT | (wavefunction[qsubset[i + 1]], wavefunction[qsubset[0]])


def inversed_CNOT_entangler(wavefunction, qsubset):
    CNOT | (wavefunction[qsubset[len(qsubset) - 1]], wavefunction[qsubset[0]])
    for i in reversed(range(len(qsubset) - 1)):
        CNOT | (wavefunction[qsubset[i]], wavefunction[qsubset[i + 1]])


def XY_full_rotation(wavefunction, qsubset, parameter):
    """
    Apply Rx(t1)Ry(t2) rotation to all the qubits with different angle
    """
    n_qubit = len(qsubset)
    for i in range(len(qsubset)):
        Rx(parameter[i]) | wavefunction[qsubset[i]]
        Ry(parameter[n_qubit + i]) | wavefunction[qsubset[i]]


def inversed_XY_full_rotation(wavefunction, qsubset, parameter):
    """
    The inversed operation of XY_full_rotation() of the same parameter
    """
    n_qubit = len(qsubset)
    for i in range(len(qsubset)):
        Ry(-parameter[n_qubit + i]) | wavefunction[qsubset[i]]
        Rx(-parameter[i]) | wavefunction[qsubset[i]]


def full_rotation(wavefunction, qsubset, parameter):
    """
    Apply Rx(t1)Ry(t2)Rx(t3) rotation to all the qubits with different angle
    """
    n_qubit = len(qsubset)
    for i in range(len(qsubset)):
        Rx(parameter[i]) | wavefunction[qsubset[i]]
        Rz(parameter[n_qubit + i]) | wavefunction[qsubset[i]]
        Rx(parameter[2 * n_qubit + i]) | wavefunction[qsubset[i]]


def inversed_full_rotation(wavefunction, qsubset, parameter):
    """
    The inversed operation of full_rotation() of the same parameter
    """
    n_qubit = len(qsubset)
    for i in range(len(qsubset)):
        Rx(-parameter[2 * n_qubit + i]) | wavefunction[qsubset[i]]
        Rz(-parameter[n_qubit + i]) | wavefunction[qsubset[i]]
        Rx(-parameter[i]) | wavefunction[qsubset[i]]


def apply_global_phase(wavefunction, angle):
    from projectq.ops import Ph
    Ph(angle) | wavefunction[0]


def generalized_rotation(wavefunction, qsubset, pauliword, evolution_time):
    """Apply e^{iPt} on the wavefunction
    Args:
        qsubset: The subset of qubits in the wavefunction that the operation applies on
        pauliword: The Pauli word P in e^{iPt} 1:X,2:Y,3:Z
        evolution_time: t in e^{iPt}
    """
    HALF_PI = math.pi / 2
    n_qubit = len(qsubset)

    # Exponentiating each Pauli string requires five parts

    # 1. Perform basis rotations
    for p in range(0, n_qubit):
        pop = pauliword[p]  # Pauli op
        if pop == 1:
            H | wavefunction[qsubset[p]]  # Hadamard
        elif pop == 2:
            Rx(HALF_PI) | wavefunction[qsubset[p]]

    # 2. First set CNOTs
    prev_index = None

    for p in range(0, n_qubit):
        pop = pauliword[p]  # Pauli op
        if pop == 0:
            continue
        if prev_index is not None:
            CNOT | (wavefunction[prev_index], wavefunction[qsubset[p]])
        prev_index = qsubset[p]

    # 3. Rotation (Note kexp & Ntrot)
    Rz(evolution_time * 2) | wavefunction[prev_index]

    # 4. Second set of CNOTs
    prev_index = None
    for p in reversed(range(0, n_qubit)):
        pop = pauliword[p]  # Pauli op
        if pop == 0:
            continue
        if prev_index is not None:
            CNOT | (wavefunction[qsubset[p]], wavefunction[prev_index])
        prev_index = qsubset[p]

    # 5. Rotate back to Z basis
    for p in range(0, len(pauliword)):
        pop = pauliword[p]  # Pauli op
        if pop == 1:
            H | wavefunction[qsubset[p]]  # Hadamard
        elif pop == 2:
            Rx(-HALF_PI) | wavefunction[qsubset[p]]
