import networkx as nx
from openfermion.ops import QubitOperator
import matplotlib.pyplot as plt


def random_graph(d=3, n=8, seed=None):
    '''
    args:   
        d (int) – The degree of each node.
        n (integer) – The number of nodes. The value of n×d must be even.
        seed (integer, random_state, or None (default)) – 
            Indicator of random number generation state.
    return: random networkx regular graph with weights. Will also print the graph.
    '''
    graph = nx.random_regular_graph(d, n, seed)
    graph.edges.data("weight", default=1)
    for _u, _v, d in graph.edges(data=True):
        d['weight'] = 1.0
    return graph


def get_maxcut_hamiltonian_from_graph(graph):
    hamiltonian = QubitOperator()
    for (i, j, wt) in graph.edges.data('weight'):
        hamiltonian += wt * QubitOperator("Z" + str(i) + " Z" + str(j))

    return hamiltonian


def get_tsp_hamiltonian_from_graph(graph):
    # Here nodes coresponding to cities, s coresponding to step
    hamiltonian = QubitOperator()
    nodes = len(graph.nodes)
    for s in range(nodes):
        for (i, j, wt) in graph.edges.data('weight'):
            hamiltonian += -wt * QubitOperator("Z" + str(i * nodes + s))
            hamiltonian += -wt * QubitOperator("Z" + str(j * nodes + s + 1))
            hamiltonian += wt * QubitOperator("Z" + str(i * nodes + s) + " Z" + str(j * nodes + s + 1))
    return hamiltonian


def get_random_maxcut_hamiltonian(d=3, n=8, seed=None):
    g = random_graph(d, n, seed)
    hamiltonian = get_maxcut_hamiltonian_from_graph(g)
    return hamiltonian


def get_random_tsp_hamiltonian(d=3, n=8, seed=None):
    g = random_graph(d, n, seed)
    hamiltonian = get_tsp_hamiltonian_from_graph(g)
    return hamiltonian
