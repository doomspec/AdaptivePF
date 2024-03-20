from ._objective import Objective, CostFunction
from ..block import HartreeFockInitBlock
from ..block._utilities import get_circuit_energy, get_state_energy
from ..backend import QuantumState, get_backend_operator


class EnergyObjective(Objective):
    """
    Attributes:
        Hamiltonian
        init_block: the block usually used for initiate the wave fucntion for the Hamiltonian.
        A block that produce a state near the ground state should be adopted.
        Usually, for the molecule, the Hartree Fock initialization is used
        obj_info: The dict for additional information of the Hamiltonian, e.g. the HF energy and the ground state energy
    Methods:
        get_cost(): generate a EnergyCost object which can be used by optimizers
    """

    def __init__(self, hamiltonian, n_qubit, init_block=None, obj_info={}):
        self.hamiltonian = hamiltonian
        self.n_qubit = n_qubit
        self.obj_info =  {"start_cost": None, "terminate_cost": None, "system_name": "Untitled"}
        self.obj_info.update(obj_info)
        self.init_block = None
        if init_block is not None:
            self.init_block = init_block
        else:
            self.init_block = HartreeFockInitBlock([])
        return

    def get_cost(self):
        return EnergyCost(self.hamiltonian)


class EnergyCost(CostFunction):

    def __init__(self, hamiltonian):
        self.hamiltonian = hamiltonian

    def get_cost_obj(self, circuit):

        ansatz = circuit.get_ansatz_on_active_parameters()
        backend_hamil = get_backend_operator(self.hamiltonian)

        def obj(parameter):
            state=ansatz(parameter)
            e=backend_hamil.get_expectation_value(state).real
            return e

        return obj

    def get_cost_value(self, circuit):
        return get_circuit_energy(circuit, self.hamiltonian)

    def get_cost_value_by_state(self, state: QuantumState):
        return get_state_energy(state, self.hamiltonian)
