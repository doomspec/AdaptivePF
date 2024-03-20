from . import inner_product


def calc_fidelity(state1,state2):
    inp=inner_product(state1,state2)
    return abs(inp)**2
def calc_infidelity(state1,state2):
    inp=inner_product(state1,state2)
    return 1-abs(inp)**2



