from projectq import MainEngine
from projectq.ops import X, All, Measure, CNOT, Z, Rz, Rx, C, TimeEvolution, Ph
from projectq.backends import CommandPrinter, CircuitDrawer
from projectq.backends import Simulator as projectq_simulator
from projectq.meta import Loop, Compute, Uncompute, Control

from projectq.cengines import (MainEngine,
                               AutoReplacer,
                               LocalOptimizer,
                               TagRemover,
                               DecompositionRuleSet)
import projectq.setups.decompositions

from time import time

"""
This file provides functions for common tasks that use the quantum simulation backend, including expectation value and amplitude
The functions accept a function called *ansatz(parameter,wavefunction)* as input.
"""


# This part can easily change to use HiQ.
# We use projectq here because it can be easily installed by pip and thus easily for Mizore to be installed.


def get_quantum_engine():
    # Create a main compiler engine with a simulator backend:
    backend = projectq_simulator(gate_fusion=True)
    cache_depth = 10
    rule_set = DecompositionRuleSet(modules=[projectq.setups.decompositions])
    engines = [TagRemover(),
               LocalOptimizer(cache_depth),
               AutoReplacer(rule_set)]
    compiler_engine = MainEngine(backend=backend, engine_list=engines)

    return compiler_engine


def evaluate_ansatz_expectation(parameter, n_qubit, hamiltonian, ansatz):
    """
    Args:

    Returns: 

    """
    from qulacs import QuantumCircuit, QuantumState
    from qulacs.quantum_operator import create_quantum_operator_from_openfermion_text

    # Initialize the wavefunction
    circuit = QuantumCircuit(n_qubit)

    # Apply the circuit
    ansatz(parameter, circuit)

    # Use the engine to implement the gates
    state = QuantumState(n_qubit)
    circuit.update_quantum_state(state)
    qulacs_hamiltonian = create_quantum_operator_from_openfermion_text(str(hamiltonian))
    return qulacs_hamiltonian.get_expectation_value(state).real


def number2bitstring(n_qubit, number):
    bitstring = [False] * n_qubit
    for i in range(0, n_qubit):
        if int((number % (1 << (i + 1))) / (1 << i)) == 1:
            bitstring[i] = True
    return bitstring


def get_ansatz_complete_amplitudes(n_qubit, ansatz):
    compiler_engine = get_quantum_engine()
    wavefunction = compiler_engine.allocate_qureg(n_qubit)
    ansatz([0] * 100, wavefunction)
    compiler_engine.flush()
    n_amp = 2 ** n_qubit
    raw_amps = [0] * n_amp
    for i in range(0, n_amp):
        raw_amps[i] = compiler_engine.backend.get_amplitude(
            number2bitstring(n_qubit, i), wavefunction)
    All(Measure) | wavefunction
    compiler_engine.flush()
    return raw_amps


def evaluate_ansatz_amplitudes(n_qubit, ansatz, bit_string_list, parameter=None):
    """
    Return the amplitudes of the kets in the bit_string_list
    Each bit_string should be a list of booleans like [False]*n_qubit
    False: |0>; True: |1>
    """
    compiler_engine = get_quantum_engine()

    # Initialize the wavefunction
    wavefunction = compiler_engine.allocate_qureg(n_qubit)

    # Apply the circuit
    if parameter is None:
        ansatz([0] * 150, wavefunction)
    else:
        ansatz(parameter, wavefunction)

    # Use the engine to implement the gates
    compiler_engine.flush()

    amp_list = []
    # Evaluate the amplitude
    for bit_string in bit_string_list:
        amp = compiler_engine.backend.get_amplitude(bit_string, wavefunction)
        amp_list.append(amp)

    # For deallocate qubit
    All(Measure) | wavefunction
    compiler_engine.flush()

    return amp_list


def evaluate_ansatz_0000_amplitudes(n_qubit, ansatz, parameter=None):
    return evaluate_ansatz_amplitudes(n_qubit, ansatz, [[False] * n_qubit], parameter)[0]


def evaluate_ansatz_1DMs(parameter, n_qubit, ansatz):
    from .WaveLocalProperties import get_one_DMs
    import numpy as np
    from openfermion.ops import QubitOperator
    compiler_engine = get_quantum_engine()
    wavefunction = compiler_engine.allocate_qureg(n_qubit)
    ansatz(parameter, wavefunction)
    compiler_engine.flush()
    one_DMs = get_one_DMs(
        compiler_engine.backend.get_expectation_value, wavefunction)

    All(Measure) | wavefunction
    compiler_engine.flush()
    return one_DMs


def evaluate_ansatz_2DMs(parameter, n_qubit, ansatz):
    from ..utilities.WaveLocalProperties import get_two_DMs
    import numpy as np
    from openfermion.ops import QubitOperator
    compiler_engine = get_quantum_engine()
    wavefunction = compiler_engine.allocate_qureg(n_qubit)
    ansatz(parameter, wavefunction)
    compiler_engine.flush()
    two_DMs = get_two_DMs(
        compiler_engine.backend.get_expectation_value, wavefunction)

    All(Measure) | wavefunction
    compiler_engine.flush()
    return two_DMs
