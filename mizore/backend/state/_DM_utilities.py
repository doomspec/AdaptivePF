from scipy.linalg import sqrtm
from numpy import matmul,trace
from math import sqrt
from .. import DensityMatrix

def calc_DM_fidelity(dm1:DensityMatrix,dm2:DensityMatrix):
    mat1=dm1.get_matrix()
    mat2=dm2.get_matrix()
    mat1=sqrtm(mat1)
    mat3=matmul(mat1,mat2)
    mat3=matmul(mat3,mat1)
    mat3=sqrtm(mat3)
    fi=abs(trace(mat3))
    fi=fi**2
    return fi