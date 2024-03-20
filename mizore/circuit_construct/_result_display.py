import os, pickle, json
from ..utilities.Visualization import draw_x_y_line_relation

RESULT_PATH = "mizore_results/"
CIRCUIT_PATH = RESULT_PATH + "block_circuits/"
CONSTRUCTION_PATH = RESULT_PATH + "contruction_runs/"
from ._circuit_constructor import NOT_DEFINED, CircuitConstructor


def mkdir(path):
    is_dir_exists = os.path.exists(path)
    if not is_dir_exists:
        os.makedirs(path)
        return True
    else:
        return False


def save_circuit(circuit, circuit_name):
    mkdir(CIRCUIT_PATH)
    path = CIRCUIT_PATH + circuit_name + ".bc"
    with open(path, "wb") as f:
        pickle.dump(circuit, f)


def save_construction(constructor, project_name):
    """
    Generate a log of circuit construction as a JSON file
    """
    save_path = CONSTRUCTION_PATH + project_name + '/'
    mkdir(save_path)
    log_dict = {}
    log_dict["project_name"] = constructor.project_name
    log_dict["pool_size"] = len(constructor.block_pool)
    log_dict["cost_list"] = constructor.cost_list
    log_dict["run_time_list"] = constructor.run_time_list
    log_dict["important_log_list"] = constructor.important_log_list
    path = save_path + "run_log.json"
    with open(path, "w") as f:
        json.dump(log_dict, f)
    path = save_path + "circuit.bc"
    with open(path, "wb") as f:
        pickle.dump(constructor.circuit, f)
    draw_current_result(constructor, save_path)


def draw_current_result(constructor: CircuitConstructor, path):
    """
    Draw the result of constrction. Two figures will be generated:
    1. Cost vs. Iteration
    2. Time vs. Iteration
    """
    mkdir(path)
    start_cost = max(constructor.cost_list)
    end_cost = 0
    red_line = None

    if constructor.terminate_cost != NOT_DEFINED:
        end_cost = constructor.terminate_cost
        red_line = end_cost
    else:
        end_cost = min(constructor.cost_list)

    space = (start_cost - end_cost) * 0.05

    draw_x_y_line_relation(range(0, len(constructor.cost_list)), constructor.cost_list, x_label="Number of Iterations",
                           y_label="Cost", filename=path + 'cost', ylim=(end_cost - space, start_cost + space),
                           red_line=red_line)

    draw_x_y_line_relation(range(0, len(constructor.run_time_list)), constructor.run_time_list,
                           x_label="Number of Iterations",
                           y_label="Time used", filename=path + 'time')
