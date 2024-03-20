from openfermion import QubitOperator

def transverse_field_ising(n_col, n_row, field=1, spin_coupling=0, periodic=True):

    n_qubit = n_col * n_row
    spin_coupling_hamiltonian = QubitOperator()

    for i in range(n_row):
        for j in range(n_col - 1):
            temp = i * n_col + j
            spin_coupling_hamiltonian += QubitOperator("Z" + str(temp)) * QubitOperator("Z" + str(temp + 1))

    for i in range(n_row - 1):
        for j in range(n_col):
            spin_coupling_hamiltonian += QubitOperator("Z" + str(i * n_col + j)) * QubitOperator(
                "Z" + str((i + 1) * n_col + j))

    if periodic:
        for i in range(n_row):
            spin_coupling_hamiltonian += QubitOperator("Z" + str(i * n_col)) * QubitOperator(
                "Z" + str((i + 1) * n_col - 1))
        temp = (n_row - 1) * n_col
        for i in range(n_col):
            spin_coupling_hamiltonian += QubitOperator("Z" + str(temp + i)) * QubitOperator("Z" + str(i))

    field_hamitonian = QubitOperator()
    for i in range(n_qubit):
        field_hamitonian += QubitOperator("X" + str(i))

    hamiltonian = spin_coupling * spin_coupling_hamiltonian + field * field_hamitonian

    return hamiltonian