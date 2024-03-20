from openfermion import QubitOperator

def get_operator_nontrivial_term_weight_sum(operator: QubitOperator):
    coeff_sum=0
    for pauli_and_coeff in operator.get_operators():
        for string_pauli in pauli_and_coeff.terms:
            if len(string_pauli) == 0:
                print("const term",string_pauli)
                continue
            coeff_sum+=abs(pauli_and_coeff.terms[string_pauli])
    return coeff_sum