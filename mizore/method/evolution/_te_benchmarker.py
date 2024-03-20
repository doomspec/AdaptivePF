from mizore.backend.state import calc_infidelity
from mizore.block import TimeEvolutionBlock,TrotterEvolutionBlock, Block
import pickle, os, json
from pathlib import Path
from mizore.backend.state import save_state , read_state

class TimeEvolutionBenchmarker():

    def __init__(self,circuits,delta_t,hamiltonian,save_path):
        self.circuits=circuits
        self.delta_t=delta_t
        self.hamiltonian=hamiltonian
        self.save_path= save_path+"/benchmark"
        Path(self.save_path).mkdir(parents=False, exist_ok=True)
        self.n_qubit=circuits[0].n_qubit
        self.n_circuit=len(self.circuits)
        self.init_circuit=circuits[0]

        with open(self.save_path + "/infidelity.json", "w") as f:
            json.dump(dict(),f)
        Path(self.save_path+"/main").mkdir(exist_ok=True)
        for i in range(0, self.n_circuit):
            save_state(self.circuits[i].get_quantum_state(), self.save_path+"/main", i)

    def generate_init_states(self):
        save_path = self.save_path + "/init"
        self.generate_benchmark_state(Block(), save_path)

    def generate_trotter_states(self,trotter_setting):
        n_trotter_step=trotter_setting["n_trotter_step"]
        n_delta_t=trotter_setting["n_delta_t"]
        trotter_block=TrotterEvolutionBlock(self.hamiltonian,n_trotter_step)([self.delta_t*n_delta_t])
        save_path=self.save_path + "/trotter{}-{}".format(n_delta_t,n_trotter_step)
        self.generate_benchmark_state(trotter_block,save_path,n_delta_t=n_delta_t)

    def generate_exact_states(self):
        te_block=TimeEvolutionBlock(self.hamiltonian)([self.delta_t])
        save_path = self.save_path + "/exact"
        self.generate_benchmark_state(te_block,save_path)

    def generate_benchmark_state(self,evolver,save_path,n_delta_t=1):
        Path(save_path).mkdir(exist_ok=True)
        generate_benchmark_state(self.init_circuit.get_quantum_state(), self.n_circuit, evolver,save_path,n_delta_t=n_delta_t)

    def calc_infidelity_list(self):
        is_exact_exist = check_complete(self.save_path,"exact",self.n_circuit)
        assert is_exact_exist
        with open(self.save_path + "/infidelity.json","r") as f:
            benchmark_dict=json.load(f)
        exist_keys=benchmark_dict.keys()
        new_keys=[]
        for folder in Path(self.save_path).iterdir():
            if not folder.is_dir():
                continue
            folder_name=folder.name
            if folder_name == "exact":
                continue
            if folder_name not in exist_keys:
                new_keys.append(folder_name)
            else:
                is_complete=check_complete(self.save_path,folder_name,self.n_circuit)
                assert is_complete

        for key in new_keys:
            infidelity_list=calc_infidelity_list(self.n_circuit,self.save_path,key)
            benchmark_dict[key]=infidelity_list

        with open(self.save_path + "/infidelity.json","w") as f:
            json.dump(benchmark_dict,f)

    def benchmark(self,trotter_settings=({"n_delta_t":10,"n_trotter_step":5},)):
        self.generate_exact_states()
        for trotter_setting in trotter_settings:
            self.generate_trotter_states(trotter_setting)
        self.generate_init_states()
        self.calc_infidelity_list()

def get_project_quantum_state_dict(save_path,n_delta_t=1,key2run=tuple())->dict:
    state_dict={}
    for folder in Path(save_path+"/benchmark/").iterdir():
        if not folder.is_dir():
            continue
        run_name = folder.name
        if run_name in key2run:
            state_dict[run_name]=get_folder_quantum_state_list(save_path,run_name,n_delta_t)
    return state_dict

def get_folder_quantum_state_list(save_path, run_name, n_delta_t=1):
    i=0
    state_list=[]
    folder_path="{}/benchmark/{}/".format(save_path,run_name)
    while True:
        if not (Path(folder_path)/"{}.wfn".format(i)).is_file():
            break
        state_list.append(read_state(folder_path,i))
        #print("{}{}.wfn".format(folder_path,i))
        i+=n_delta_t
    return state_list


def check_complete(save_path,benchmark_name,n_circuit):
    return (Path(save_path)/ benchmark_name / "{}.wfn".format(n_circuit-1)).is_file()


def calc_infidelity_list(n_circuit,save_path,benchmark_name):
    infidelity_list = [-1] * n_circuit
    for i in range(0,n_circuit):
        if not (Path(save_path) / benchmark_name / "{}.wfn".format(i)).is_file():
            continue
        circuit_state = read_state("/".join([save_path,"exact"]),i)
        ref_state = read_state("/".join([save_path,benchmark_name]),i)
        infidelity_list[i]=calc_infidelity(circuit_state,ref_state)
    return infidelity_list


def generate_benchmark_state(init_state,n_circuit,evolver,save_path,n_delta_t=1):
    #print(save_path)
    ref_state = init_state
    save_state(ref_state, save_path, 0)
    for i in range(1, n_circuit):
        if i%n_delta_t!=0:
            continue
        evolver | ref_state
        save_state(ref_state,save_path,i)
        #print(ref_state)


def read_project_for_benchmark(path,delta_n=1):
    circuits=read_circuits_from_file(path+"/circuits",delta_n=delta_n)
    with open(path+"/hamiltonian.pickle", "rb") as f:
        hamiltonian = pickle.load(f)
    with open(path+"/para_dict.json", "r") as f:
        delta_t = json.load(f)["delta_t"]
    return circuits, delta_t, hamiltonian, path


def read_circuits_from_file(path, delta_n=1):
    circuit_list = []
    circuit_index = 0
    while True:
        circuit_path = path+"/"+str(circuit_index)+".bc"
        if not os.path.exists(circuit_path):
            break
        with open(circuit_path, "rb") as f:
            step_circuit = pickle.load(f)
        circuit_list.append(step_circuit)
        circuit_index += delta_n
    return circuit_list


def benchmarker_from_path(project_path):
    project = read_project_for_benchmark(project_path)
    return TimeEvolutionBenchmarker(*project)
