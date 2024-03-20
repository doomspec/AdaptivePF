import itertools

"""
Common iterators that are used in Mizore
"""


def iter_qsubset_pauli_of_operator(operator):
    for pauli_and_coeff in operator.get_operators():
        for string_pauli in pauli_and_coeff.terms:
            if string_pauli != ():
                yield string_pauli2qsubset_pauli(string_pauli)


def iter_coeff_qsubset_pauli_of_operator(operator):
    for pauli_and_coeff in operator.get_operators():
        for string_pauli in pauli_and_coeff.terms:
            if string_pauli != ():
                coeff = pauli_and_coeff.terms[string_pauli]
                res = [coeff]
                res.extend(string_pauli2qsubset_pauli(string_pauli))
                yield res


def iter_terms_in_fermion_operator(operator):
    from openfermion.ops import FermionOperator
    for pauli_and_coeff in operator.get_operators():
        for string_pauli in pauli_and_coeff.terms:
            yield FermionOperator(string_pauli)


def iter_terms_in_qubit_operator(operator):
    from openfermion.ops import QubitOperator
    for pauli_and_coeff in operator.get_operators():
        for string_pauli in pauli_and_coeff.terms:
            yield QubitOperator(string_pauli)


def iter_partial_operators(hamiltonian, max_n_term_in_new_operators):
    from openfermion.ops import QubitOperator
    new_operator = QubitOperator()
    n_terms_added = 0
    for string_pauli in hamiltonian.terms:
        new_operator += QubitOperator(string_pauli) * hamiltonian.terms[string_pauli]
        n_terms_added += 1
        if n_terms_added == max_n_term_in_new_operators:
            yield new_operator
            new_operator = QubitOperator()
            n_terms_added = 0
    if n_terms_added != 0:
        yield new_operator


def string_pauli2qsubset_pauli(string_pauli, make_imaginary=False):
    qsubset = []
    pauli = []
    pauli_dict = {'X': 1, 'Y': 2, 'Z': 3}
    for term in string_pauli:
        qsubset.append(term[0])
        pauli.append(pauli_dict[term[1]])
    return qsubset, pauli


def is_contain_odd_Y(pauli):
    y_num = 0
    for i in pauli:
        if i == 2:
            y_num += 1
    is_valid = (y_num % 2 == 1)
    return is_valid


def iter_odd_Y_pauli(length):
    for pauli in itertools.product(range(1, 4), repeat=length):
        if not is_contain_odd_Y(pauli):
            continue
        yield pauli


def iter_qsubset(length, indices):
    for qs in itertools.combinations(range(0, len(indices)), length):
        qsubset = []
        for i in qs:
            qsubset.append(indices[i])
        yield qsubset


def iter_qsubset_odd_Y_pauli_by_length(length, indices):
    for qsubset in iter_qsubset(length, indices):
        for pauli in iter_odd_Y_pauli(length):
            yield qsubset, pauli


def iter_qsubset_odd_Y_pauli(indices):
    for length in range(len(indices)):
        for qsubset, pauli in iter_qsubset_odd_Y_pauli_by_length(length + 1, indices):
            yield qsubset, pauli


def iter_all_qsubset_pauli_by_length(length, indices):
    for qsubset in iter_qsubset(length, indices):
        for pauli in itertools.product(range(1, 4), repeat=length):
            yield qsubset, pauli


if __name__ == "__main__":
    from openfermion.ops import FermionOperator

    a = FermionOperator('4^ 3 9 3^') + FermionOperator('4^  3^')
    for term in iter_terms_in_fermion_operator(a):
        print(term)
