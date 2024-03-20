from openfermion.transforms import bravyi_kitaev
from openfermion.transforms import jordan_wigner
from openfermion.transforms import binary_code_transform, parity_code
from openfermion.ops import FermionOperator, QubitOperator

"""
Methods for generating usually used fermion-qubit transformation
"""


def get_parity_transform(n_spinorbital):
    def parity_transform(FermionOperator):
        return binary_code_transform(FermionOperator, parity_code(n_spinorbital))

    return parity_transform


def make_transform_spin_separating(transform, n_spinorbital):
    def new_transform(fermion_operator: FermionOperator):
        separated_operator = separate_odd_even(fermion_operator, n_spinorbital)
        return transform(separated_operator)

    return new_transform


def separate_odd_even(fermion_operator, n_spinorbital):
    return rearrange_sites(fermion_operator, separate_odd_even_map(n_spinorbital))


def rearrange_sites(fermion_operator, change_rule):
    """
    Can only deal with fermion and qubit operator!!
    Rearrange the order of the spin orbitals
    -- change_rule: should be an array e.g. [0,3,2,1,4,5,6,7]  means exchange 3,1
    """
    new_fermion_operator = None
    is_qubit_operator = True
    if isinstance(fermion_operator, QubitOperator):
        new_fermion_operator = QubitOperator()
        is_qubit_operator = True
    if isinstance(fermion_operator, FermionOperator):
        new_fermion_operator = FermionOperator()
        is_qubit_operator = False

    for pauli_and_coff in fermion_operator.get_operators():
        for string_pauli in pauli_and_coff.terms:
            # string_pauli is like ( (1, 'X'), (2, 'Y')  )
            new_str = list(string_pauli)
            for i in range(0, len(string_pauli)):
                new_str[i] = (change_rule[string_pauli[i][0]], string_pauli[i][1])

            if is_qubit_operator:
                new_fermion_operator += pauli_and_coff.terms[string_pauli] * \
                                        QubitOperator(new_str)
            else:
                new_fermion_operator += pauli_and_coff.terms[string_pauli] * \
                                        FermionOperator(new_str)

    if new_fermion_operator == None:
        raise Exception("Error, please check the type of the input fermion_operator!")

    return new_fermion_operator


def separate_odd_even_map(n_spinorbital):
    """
    Put the spin orbital of odd index to the last, making the first half spin up and follows spin down
    """
    if n_spinorbital % 2 != 0:
        print("n_spinorbitals is not even!")
        return
    rule = [0] * n_spinorbital
    for i in range(0, n_spinorbital):
        if i % 2 == 0:
            rule[i] = i // 2
        else:
            rule[i] = i // 2 + n_spinorbital // 2
    return rule


if __name__ == "__main__":
    a = FermionOperator("0^ 1^ 2 3")
    parity = get_parity_transform(4)
    print(parity(a))
    new_transform = make_transform_spin_separating(parity, 4)
    print()
    print(new_transform(a))
