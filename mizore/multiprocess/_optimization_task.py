from ._task import Task
from ..block import BlockCircuit
from ..optimizer import ParameterOptimizer
from openfermion.ops import QubitOperator
from ..objective._objective import CostFunction


class OptimizationTask(Task):
    """
    The task for finding the optimized cost and parameter for a Block Circuit
    The result is usually a tuple (cost,amp)
    See the document of optimizers to see what exact the result is
    """

    def __init__(self, circuit: BlockCircuit, optimizer: ParameterOptimizer, cost: CostFunction):
        Task.__init__(self)
        self.circuit = circuit
        self.optimizer = optimizer
        self.cost = cost

    def run(self):
        res = self.optimizer.run_optimization(self.circuit, self.cost)
        return res
