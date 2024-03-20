from ._trotter_block import TrotterEvolutionBlock
from ..backend.gate import PauliRotation, Probabilistic, DenseMatrix
from openfermion.ops import QubitOperator
from numpy import sign
import numpy as np
from numpy.linalg import matrix_power
from ..utilities.Tools import get_operator_qsubset

class QDriftBlock(TrotterEvolutionBlock):

    def __init__(self, hamiltonian: QubitOperator, n_trotter_step=1, evolution_time=0):
        TrotterEvolutionBlock.__init__(self,hamiltonian,n_trotter_step=n_trotter_step,evolution_time=evolution_time)
        #self.qsubset = get_operator_qsubset(hamiltonian)
        self.n_term=len(self.qsubset_pauliword_list)
        coeff_list = np.array([self.qsubset_pauliword_list[i][0] for i in range(0,self.n_term)])
        self.probability_list = abs(coeff_list)
        self.sign_list = sign(coeff_list)
        self.interaction_strength= sum(self.probability_list)
        self.probability_list = self.probability_list/self.interaction_strength

    def get_forward_gates(self, parameter):
        gate_list = [0]*self.n_term
        temp=-2*self.interaction_strength*(parameter[0]+self.parameter[0])/self.n_trotter_step
        for i in range(0,self.n_term):
            gate_list[i]=PauliRotation(self.qsubset_pauliword_list[i][1],self.qsubset_pauliword_list[i][2],
                                       temp*self.sign_list[i])
            #print(self.qsubset_pauliword_list[i][1],self.qsubset_pauliword_list[i][2],self.qsubset_pauliword_list[i][0],temp*self.sign_list[i])
        prob_gate=Probabilistic(self.probability_list,gate_list)
        return [prob_gate]*self.n_trotter_step

    def get_inverse_gates(self, parameter):
        assert False
        pass

    def get_gate_used(self):
        gate_used=TrotterEvolutionBlock.get_gate_used(self)
        for key,value in gate_used.items():
            gate_used[key]=value/self.n_term
        return gate_used

    def __str__(self):
        info = self.basic_info_string()
        info += "; qDRIFT_TimeEvolution: T=" + \
                str(self.parameter[0]) + " N_step:" + str(self.n_trotter_step)
        return info