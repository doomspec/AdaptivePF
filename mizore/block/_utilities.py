from ._block_circuit import BlockCircuit
from copy import copy


def concatenate_circuit(first_circuit: BlockCircuit = None, second_circuit: BlockCircuit = None):
    return concatenate_circuit_list((first_circuit, second_circuit))


def concatenate_circuit_list(circuit_list=None):
    n_qubit = circuit_list[0].n_qubit
    new_active_qubit_list = []
    new_circuit = BlockCircuit(n_qubit)
    for circuit in circuit_list:
        if circuit.n_qubit != n_qubit:
            raise Exception(
                "The number of qubit of the circuits is not unified!")
        n_block_added = len(new_circuit.block_list)
        for active in circuit.active_position_list:
            new_active_qubit_list.append(active + n_block_added)
        for block in circuit.block_list:
            new_circuit.add_copied_block(copy(block))
        new_circuit.active_position_list = copy(new_active_qubit_list)
    return new_circuit


def get_inverse_circuit(circuit: BlockCircuit):
    new_circuit = BlockCircuit(circuit.n_qubit)
    for block in reversed(circuit.block_list):
        reversed_block = copy(block)
        reversed_block.is_inverse = True
        new_circuit.add_block_without_copy(reversed_block)
    return new_circuit


from qulacs import QuantumCircuit, QuantumState
from qulacs.quantum_operator import create_quantum_operator_from_openfermion_text


def get_circuit_energy(block_circuit: BlockCircuit, hamiltonian):
    qulacs_hamiltonian = create_quantum_operator_from_openfermion_text(str(hamiltonian))
    circuit = block_circuit.get_circuit_by_parameters([], [])
    state = QuantumState(block_circuit.n_qubit)
    circuit.update_quantum_state(state)
    return qulacs_hamiltonian.get_expectation_value(state).real


def get_state_energy(state: QuantumState, hamiltonian):
    qulacs_hamiltonian = create_quantum_operator_from_openfermion_text(str(hamiltonian))
    return qulacs_hamiltonian.get_expectation_value(state).real
