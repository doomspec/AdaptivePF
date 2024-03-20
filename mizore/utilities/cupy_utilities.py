from numpy import real,imag
from cupy.linalg import lstsq as culstsq
from cupy.core import array as cuarray

def calc_derivative_cupy(mat_A,vec_C):
    return culstsq(cuarray(real(mat_A)), cuarray(imag(vec_C)),rcond=None)[0].get()