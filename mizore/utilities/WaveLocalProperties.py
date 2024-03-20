from openfermion.ops import QubitOperator
import numpy as np
import math

PauliI_half = np.array([[0.5, 0], [0, 0.5]], np.complex)
PauliI = np.array([[1, 0], [0, 1]], np.complex)
PauliX = np.array([[0, 1], [1, 0]], np.complex)
PauliY = np.array([[0, -1j], [1j, 0]], np.complex)
PauliZ = np.array([[1, 0], [0, -1]], np.complex)
PauliMat = [PauliI, PauliX, PauliY, PauliZ]
PauliString = ['I', 'X', 'Y', 'Z']
import scipy.linalg as la

"""
The functions for analyzing the wavefunction obtained from classical pre-calculations
"""


########################
#   Density Matrices
########################

def get_one_DMs(get_expectation_value, wavefunction):
    n_qubit = len(wavefunction)
    one_DMs = np.array([PauliI] * n_qubit)
    for i in range(0, n_qubit):
        for ip in range(1, 4):
            one_DMs[i] = one_DMs[i] + get_expectation_value(
                QubitOperator(PauliString[ip] + str(0)), [wavefunction[i]]) * PauliMat[ip]
        one_DMs[i] = one_DMs[i] / 2
    return one_DMs


def get_two_DMs(get_expectation_value, wavefunction):
    n_qubit = len(wavefunction)
    two_DMs = np.array([[np.kron(PauliI, PauliI)] * n_qubit] * n_qubit)
    for i in range(0, n_qubit):
        for j in range(0, i):
            for ip in range(0, 4):
                for jp in range(0, 4):
                    if ip == 0 and jp == 0:
                        continue
                    if ip == 0:
                        itr_pauli_operator = QubitOperator(
                            PauliString[jp] + str(j))
                    if jp == 0:
                        itr_pauli_operator = QubitOperator(
                            PauliString[ip] + str(i))
                    if ip != 0 and jp != 0:
                        itr_pauli_operator = QubitOperator(
                            PauliString[ip] + str(i) + ' ' + PauliString[jp] + str(j))
                    two_DMs[i][j] = two_DMs[i][j] + get_expectation_value(
                        itr_pauli_operator, wavefunction) * np.kron(PauliMat[ip], PauliMat[jp])
            two_DMs[i][j] = two_DMs[i][j] / 4
    return two_DMs


def two_DM_to_one_DMs(two_DM):
    one_DM_1 = np.array([[0.0 + 0.0j] * 2] * 2)
    one_DM_2 = np.array([[0.0 + 0.0j] * 2] * 2)
    one_DM_1[0][0] = two_DM[0][0] + two_DM[1][1]
    one_DM_1[1][0] = two_DM[2][0] + two_DM[3][1]
    one_DM_1[0][1] = two_DM[0][2] + two_DM[1][3]
    one_DM_1[1][1] = two_DM[2][2] + two_DM[3][3]

    one_DM_2[0][0] = two_DM[0][0] + two_DM[2][2]
    one_DM_2[1][0] = two_DM[0][1] + two_DM[2][3]
    one_DM_2[0][1] = two_DM[1][0] + two_DM[3][2]
    one_DM_2[1][1] = two_DM[1][1] + two_DM[3][3]
    return one_DM_1, one_DM_2


def get_processed_eigv_for_DM(density_mat):
    eigv = np.linalg.eigvalsh(density_mat)
    # Round the eigenvalues
    realeigv = np.round(np.real(eigv), 7)
    # Renormalize the eigenvalues
    norm = 0
    for i in range(0, len(density_mat[0])):
        norm += realeigv[i]
    if norm == 0:
        print('Density matrix not proper, eigenvalue sum to 0, 0 returned')
        return 0
    else:
        realeigv = realeigv / norm
    return realeigv


#####################################
#   Entropy and Mutual Information
#####################################
def entropy_one_DMs(one_DMs):
    return [entropy_one_DM(one_DMs[i]) for i in range(len(one_DMs))]


def entropy_one_DM(one_DM: np.array):
    realeigv = get_processed_eigv_for_DM(one_DM)
    entropy = 0
    for i in range(0, len(one_DM[0])):
        if realeigv[i] != 0 and realeigv[i] != 1:
            assert realeigv[i] > 0
            entropy += -realeigv[i] * math.log2(realeigv[i])
    return entropy


def purity_one_DM(one_DM: np.array):
    realeigv = get_processed_eigv_for_DM(one_DM)
    purity = 0
    for i in range(0, len(one_DM[0])):
        purity += realeigv[i] * realeigv[i]
    return purity


def get_mutual_information_by_2DMs(two_DMs):
    n_qubit = len(two_DMs)
    mutual_information = np.array([[0.0] * n_qubit] * n_qubit)
    for i in range(0, n_qubit):
        for j in range(0, i):
            one_DM_1, one_DM_2 = two_DM_to_one_DMs(two_DMs[i][j])
            mutual_information[i][j] = 0.5 * (
                    entropy_one_DM(one_DM_1) + entropy_one_DM(one_DM_2) - entropy_one_DM(two_DMs[i][j]))
            mutual_information[j][i] = mutual_information[i][j]
    return mutual_information


#################################
#   Entanglement of Formation
#################################

# Concurrence
def Concurrence(rho):
    sy = np.array([[0., -1j], [1j, 0.]])
    sqrtRho = la.sqrtm(rho)
    YY = np.kron(sy, sy)
    matrixR = sqrtRho * YY * np.conjugate(rho) * YY * sqrtRho
    eigvals, eigvecs = la.eigh(matrixR)
    eigvals.sort()  # sort eigenvalues from smallest to largest
    return max(0, 4 * (eigvals[3] - eigvals[2] - eigvals[1] - eigvals[0]))


# Entanglement of formation
def EoF(rho):
    # check that the input 'rho' is of 4-by-4 matrix and Hermitian
    rows = len(rho)
    columns = len(rho[0])
    rhoC = np.conjugate(rho)
    if (rows == 4 and columns == 4):
        C = Concurrence(rho)
        plus = 0.5 * (1 + np.sqrt(1 - C ** 2)) * \
               np.lib.scimath.log2(0.5 * (1 + np.sqrt(1 - C ** 2)))
        if (0.5 * (1 - np.sqrt(1 - C ** 2)) != 0):  # to avoid undefined when C = 0
            minus = 0.5 * (1 - np.sqrt(1 - C ** 2)) * \
                    np.lib.scimath.log2(0.5 * (1 - np.sqrt(1 - C ** 2)))
        else:
            minus = 0.5 * (1 - np.sqrt(1 - C ** 2))
        return -plus - minus
    else:
        print("The input density matrix is not two-qubit type or not Hermitian or both. \
        Please enter the correct one.")


def EoF_regularize(EoF):
    if EoF <= 0:
        return 0
    else:
        EoF = 1 / (math.log2(EoF) ** 6)
        return EoF


def get_EoF_by_2DMs(twodensity):
    n_qubit = len(twodensity)
    entanglement_of_formation = np.array([[0.0] * n_qubit] * n_qubit)
    for i in range(0, n_qubit):
        for j in range(0, i):
            entanglement_of_formation[i][j] = EoF_regularize(EoF(twodensity[i][j]))
            entanglement_of_formation[j][i] = entanglement_of_formation[i][j]
    return entanglement_of_formation
