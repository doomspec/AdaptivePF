from pickle import load,dump
from .._utilities import QuantumStateFromVec
def save_state(state, save_path, file_name):
    with open(save_path + "/{}.wfn".format(file_name), "wb") as f:
        dump(state.get_vector(), f)

def read_state(save_path, file_name):
    return  QuantumStateFromVec(read_state_vec(save_path, file_name))

def read_state_vec(save_path, file_name):
    with open(save_path + "/{}.wfn".format(file_name),"rb") as f:
        state_vec=load(f)
    return  state_vec