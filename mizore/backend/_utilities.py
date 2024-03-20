from . import QuantumCircuit, QuantumState
from qulacs.quantum_operator import create_quantum_operator_from_openfermion_text
from qulacs import GeneralQuantumOperator
import numpy as np
from numpy.linalg import norm

def get_circuit_by_gate_list(n_qubit, gate_list):
    circuit = QuantumCircuit(n_qubit)
    for gate in gate_list:
        circuit.add_gate(gate)
    return circuit

def get_backend_operator(openfermion_op)->GeneralQuantumOperator:
    return create_quantum_operator_from_openfermion_text(str(openfermion_op))

from math import log2
def QuantumStateFromVec(state_vec)-> QuantumState:
    n_qubit=int(log2(len(state_vec)))
    state=QuantumState(n_qubit)
    state.load(state_vec)
    return state

def QuantumStateFromGeneralVec(state_vec)-> QuantumState:
    n_qubit=int(log2(len(state_vec)))
    state_vec=np.array(state_vec)
    vec_norm=norm(state_vec)
    state_vec = state_vec/vec_norm
    state=QuantumState(n_qubit)
    state.load(state_vec)
    return state
