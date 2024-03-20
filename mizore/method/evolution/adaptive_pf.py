from mizore.method.evolution import *
from mizore.block import BlockCircuit,GlobalPhaseBlock
from mizore.block_pool._rotation_pools import hamiltonian_rotations_pool
from mizore.block_pool import BlockPool
from mizore.multiprocess import TaskManager

def run_adaptive_product_formula(energy_obj,quality_cutoff=0.3,n_circuit=30,delta_t=0.2,stepsize=2e-3,n_worker=40,nM_cutoff=5e8,save_root="/home/mh-group/mizore_results/adaptive_evolution"):

    hamil=energy_obj.hamiltonian
    n_qubit=energy_obj.n_qubit
    bc=BlockCircuit(n_qubit)
    bc.add_copied_block(GlobalPhaseBlock())
    bc.add_copied_block(energy_obj.init_block)
    pool = BlockPool(hamiltonian_rotations_pool(hamil))
    bc.set_all_block_active()
    with TaskManager(n_worker) as tm:
        evolver=TimeEvolutionEvolver(bc.n_qubit,hamil,quality_cutoff=quality_cutoff,stepsize=stepsize,nM_cutoff=nM_cutoff)
        constructor=TimeEvolutionConstructor(n_qubit,pool,hamil,quality_cutoff/2,tm,nM_cutoff*0.9)
        smart_evolver=QualityEnsuredEvolver(bc,constructor,evolver,para_dict={"project_name":"{}_{}".format(energy_obj.obj_info["system_name"],str(quality_cutoff)),"is_to_save":True, "n_circuit":n_circuit, "save_root":save_root, "delta_t": delta_t})
        circuit_list=smart_evolver.run(construct_first=True)
    print("Final circuit:")
    print(circuit_list[-1])
    print("Result saved at ",smart_evolver.save_path)

    return smart_evolver.save_path