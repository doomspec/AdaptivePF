from openfermion import QubitOperator
from copy import copy

def get_random_coeff_operator(operator,start,end):
    import random
    new_operator=QubitOperator()
    for pauli_and_coeff in operator.get_operators():
        for string_pauli in pauli_and_coeff.terms:
            new_operator+=random.uniform(start,end)*QubitOperator(string_pauli)
    return new_operator

def normalize_operator_coeff(operator,average):
    import random
    coeff_sum=0
    for pauli_and_coeff in operator.get_operators():
        for string_pauli in pauli_and_coeff.terms:
            coeff_sum+=abs(pauli_and_coeff.terms[string_pauli])
    new_operator=copy(operator)*(average/coeff_sum)
    return new_operator

