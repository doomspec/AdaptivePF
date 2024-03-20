
from ..backend.gate import DenseMatrix
from . import Block
from numpy.linalg import inv
class MatrixBlock(Block):
    n_parameter = 0
    IS_INVERSE_DEFINED = True
    IS_LOCALIZE_AVAILABLE = False
    IS_DERIVATIVE_DEFINE = False

    def __init__(self, qsubset, matrix, init_angle=0):
        Block.__init__(self, n_parameter=0)
        self.qsubset=qsubset
        self.matrix=matrix
        self.inverse_matrix=inv(matrix)

    def get_forward_gates(self, parameter):
        return [DenseMatrix(self.qsubset,self.matrix)]
    def get_inverse_gates(self, parameter):
        return [DenseMatrix(self.qsubset,self.inverse_matrix)]

    def get_gate_used(self):
        return {}

    def __str__(self):
        info = self.basic_info_string()
        return info
