from ._task import Task
from ..block import BlockCircuit
from ..objective._objective import CostFunction


class EvaluationTask(Task):

    def __init__(self, circuit: BlockCircuit, cost: CostFunction):
        Task.__init__(self)
        self.circuit = circuit
        self.cost = cost

    def run(self):
        return self.cost.get_cost_value(self.circuit)
