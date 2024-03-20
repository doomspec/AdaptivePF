import numpy as np
from numpy import array,dot, vdot, asarray, conjugate, zeros, imag, real
from numpy.linalg import norm


def get_parameter_directions(circuit):
    parameter_directions = []
    for position in circuit.active_positions:
        active_block = circuit.block_list[position]
        block_n_para = active_block.n_parameter
        for in_block_position in range(block_n_para):
            derivative_block = active_block.get_derivative_block(in_block_position)
            original_block=circuit.block_list[position]
            circuit.block_list[position] = derivative_block
            parameter_directions.append(circuit.get_quantum_state().get_vector())
            circuit.block_list[position]=original_block
    return array(parameter_directions)

def calc_C_vec(parameter_directions, hamil_state_vec):
    return array([vdot(direction, hamil_state_vec) for direction in parameter_directions]) #TODO check the minus sign

def calc_A_mat(parameter_directions):
    dim=len(parameter_directions)
    mat_A=zeros((dim,dim),dtype=complex)
    for i in range(dim):
        for j in range(i):
            mat_A[i][j] = vdot(parameter_directions[i],parameter_directions[j])
            mat_A[j][i] = conjugate(mat_A[i][j])
    for i in range(dim):
        mat_A[i][i] = norm(parameter_directions[i])**2

    return mat_A

from numba import jit
from numpy.linalg import lstsq
@jit(nopython=True)
def calc_derivative(mat_A_,vec_C_):
    mat_A = real(mat_A_)#+np.eye(len(mat_A_))*1e-6
    vec_C = imag(vec_C_)
    eps = 1e-4
    lam = lstsq(mat_A, vec_C)[0]
    # print(lam)
    signs = np.sign(lam)
    lam = lstsq(mat_A, vec_C - eps * signs)[0]
    # print(lam)
    new_signs = np.sign(lam)
    max_iter = 5
    i = 0
    while norm(signs - new_signs, ord=1) != 0 and i < max_iter:
        signs = new_signs
        lam = lstsq(mat_A, vec_C - eps * signs)[0]
        new_signs = np.sign(lam)
        # print(lam)
        i += 1
    return lam
    #return lstsq(real(mat_A)+np.eye(len(mat_A))*1e-6,imag(vec_C), rcond=1e-2)[0]

def calc_derivative0(mat_A,vec_C):
    return lstsq(real(mat_A),imag(vec_C))[0]

def calc_derivative_quality(np_hamil,circuit):
    parameter_directions = get_parameter_directions(circuit)
    current_vec=array(circuit.get_quantum_state().get_vector())
    hamil_state_vec=asarray(dot(np_hamil,current_vec))[0]
    vec_C=calc_C_vec(parameter_directions,hamil_state_vec)
    mat_A=calc_A_mat(parameter_directions)
    derivative=calc_derivative(mat_A,vec_C)
    quality=calc_quality(derivative,parameter_directions,hamil_state_vec)
    return real(derivative), quality, mat_A, vec_C

def calc_quality(derivative,parameter_directions,hamil_state_vec):
    dim=len(hamil_state_vec)
    direction_to_move=zeros((dim,),dtype=complex)
    for i in range(len(derivative)):
        direction_to_move+=derivative[i]*parameter_directions[i]
    direction_diff=direction_to_move-(-1j)*hamil_state_vec
    return norm(direction_diff)