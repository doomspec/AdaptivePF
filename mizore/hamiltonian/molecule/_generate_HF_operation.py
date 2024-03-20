from openfermion.ops import FermionOperator, QubitOperator
from openfermion.transforms import bravyi_kitaev, jordan_wigner
from openfermion.utils import hermitian_conjugated

"""
The functions for generating the Hartree Fock initial operator based on the transformation and number of electrons
See get_HF_operator()
"""


def get_dressed_operator(dress_operator, operator2dress):
    """ Return a QubitOperator modified as U^HU
    Args:
        dress_operator: U
        operator2dress: H
    """
    new_operator = hermitian_conjugated(dress_operator) * operator2dress
    new_operator = new_operator * dress_operator
    return new_operator


def get_electron_fermion_operator(n_electron):
    """
    Return a fermion operator creates N electrons
    a0^ a1^ ... aN^
    Args:
        n_electron: N
    """
    creations = []
    for i in range(0, n_electron):
        creations.append((i, 1))
    fermi_operator = FermionOperator(creations)
    return fermi_operator


def get_0000_state_operator(qubit_operator):
    """
    Return a simplified QubitOperator U0 s.t.
    U0|0..00>=U1|0..00>
    Args:
        qubit_operator:U1
    """
    new_qubit_operator = QubitOperator()
    for pauli_and_coff in qubit_operator.get_operators():
        for string_pauli in pauli_and_coff.terms:
            new_string = ""
            new_coff = 1
            for terms in string_pauli:
                if terms[1] == 'X':
                    new_string += "X" + str(terms[0]) + ' '
                if terms[1] == 'Y':
                    new_coff *= 1j
                    new_string += "X" + str(terms[0]) + ' '
            new_qubit_operator += new_coff * pauli_and_coff.terms[string_pauli] * QubitOperator(new_string)
    new_qubit_operator.compress()
    return new_qubit_operator


def get_HF_operator(n_electron, transformation):
    """
    Return a QubitOperator that trans produce Hartree-Fock state on |0..00> state.
    """
    creations = []
    for i in range(0, n_electron):
        creations.append((i, 1))
    fermi_operator = FermionOperator(creations)
    # print(fermi_operator)
    qubit_operator = transformation(fermi_operator)
    # print(qubit_operator)
    qubit_operator_on_0000_state = get_0000_state_operator(qubit_operator)
    # print(qubit_operator_on_0000_state)
    return qubit_operator_on_0000_state


if __name__ == "__main__":
    x = 1j * QubitOperator("")
    y = QubitOperator("Y0") + 1j * QubitOperator("Z0")
    print(get_dressed_operator(y, x))
