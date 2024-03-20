from ._gate import QuantumGate
from mizore.backend import DensityMatrix,QuantumState
from mizore.backend._utilities import QuantumStateFromVec
class StatePreparationGate(QuantumGate):
    def __init__(self,state_vector):
        self.state_vector=state_vector
    def update_quantum_state(self,quantum_state):
        assert not isinstance(quantum_state,DensityMatrix)
        quantum_state.load(self.state_vector)
