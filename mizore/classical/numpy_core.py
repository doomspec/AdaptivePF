from .sparse_tools import get_sparse_operator
import scipy, numpy


def get_kth_excited_state(k, n_qubits, operator, initial_guess=None):
    """
    Adopted from openfermion.linalg.get_ground_state
    """
    sparse_operator = get_sparse_operator(operator, n_qubits)

    values, vectors = scipy.sparse.linalg.eigsh(sparse_operator,
                                                k=k + 1,
                                                v0=initial_guess,
                                                which='SA',
                                                maxiter=1e7)

    order = numpy.argsort(values)
    values = values[order]
    vectors = vectors[:, order]
    eigenvalue = values[k]
    eigenstate = vectors[:, k]
    eigenstate=eigenstate.T

    #print(numpy.sqrt((numpy.abs(eigenstate)**2).sum())) # To ensure normalized

    return eigenvalue, eigenstate

def get_eigen_value_states(rank_list, n_qubits, operator, initial_guess=None):
    """
    :param rank_list: ranks of the wanting eigenvalues and eigenvectors
    :return: two lists eigenvalue_list, eigenstate_list
    """
    max_order=max(rank_list)

    sparse_operator = get_sparse_operator(operator, n_qubits)

    values, vectors = scipy.sparse.linalg.eigsh(sparse_operator,
                                                k=max_order + 1,
                                                v0=initial_guess,
                                                which='SA',
                                                maxiter=1e7)

    eigenvalue_list=[None]*len(rank_list)
    eigenstate_list=[None]*len(rank_list)

    order = numpy.argsort(values)
    values = values[order]
    vectors = vectors[:, order]
    for i in range(len(rank_list)):
        eigenvalue = values[rank_list[i]]
        eigenstate = vectors[:, rank_list[i]]
        eigenstate=eigenstate.T
        eigenvalue_list[i]=eigenvalue
        eigenstate_list[i]=eigenstate
    #print(numpy.sqrt((numpy.abs(eigenstate)**2).sum())) # To ensure normalized
    return eigenvalue_list, eigenstate_list


def get_operator_matrix(operator, n_qubits=None, trunc=None, hbar=1):
    return get_sparse_operator(operator, n_qubits=n_qubits, trunc=trunc, hbar=hbar).todense()