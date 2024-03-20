from scipy.optimize import minimize, basinhopping
from ..utilities.Tools import random_list
from ._parameter_optimizer import ParameterOptimizer
from ..objective._objective import CostFunction


class BasinhoppingOptimizer(ParameterOptimizer):

    def __init__(self, random_initial=0.1, niter=10, temperature=0.1, stepsize=1e-5, tol=1e-5):
        ParameterOptimizer.__init__(self)

        self.random_initial = random_initial
        self.niter = niter
        self.temperature = temperature
        self.stepsize = stepsize
        self.tol = tol

    def run_optimization(self, circuit, cost: CostFunction):
        n_parameter = circuit.get_active_n_parameter()
        initial_parameter = [0.0] * n_parameter
        if self.random_initial != 0:
            initial_parameter = random_list(-self.random_initial, self.random_initial, n_parameter)

        obj = cost.get_cost_obj(circuit)

        opt_result = basinhopping(obj, initial_parameter, niter=self.niter,
                                  T=self.temperature, stepsize=self.stepsize, minimizer_kwargs={
                "tol": self.tol, "method": 'BFGS'},
                                  take_step=None, accept_test=None, callback=None,
                                  disp=False, niter_success=None)

        return opt_result.fun, opt_result.x
