from mizore.method.evolution import *
from mizore.block import BlockCircuit, GlobalPhaseBlock
from mizore.block_pool._rotation_pools import hamiltonian_rotations_pool
from mizore.block_pool import BlockPool
from mizore.hamiltonian import make_example_H6
from mizore.multiprocess import TaskManager


def evolve(quality_cutoff=0.3, n_circuit=30, delta_t=0.2):
    obj = make_example_H6()
    hamil = obj.hamiltonian
    n_qubit = obj.n_qubit
    bc = BlockCircuit(n_qubit)
    bc.add_copied_block(GlobalPhaseBlock())
    bc.add_copied_block(obj.init_block)
    print("n_qubit", n_qubit)
    pool = BlockPool(hamiltonian_rotations_pool(hamil))
    print(bc)
    bc.set_all_block_active()
    save_root = "/home/mh-group/"
    with TaskManager(40) as tm:
        evolver = TimeEvolutionEvolver(bc.n_qubit, hamil, quality_cutoff=quality_cutoff, stepsize=2e-3)
        constructor = TimeEvolutionConstructor(n_qubit, pool, hamil, quality_cutoff / 2, tm)
        smart_evolver = QualityEnsuredEvolver(bc, constructor, evolver,
                                              para_dict={"project_name": "H2O_" + str(quality_cutoff),
                                                         "is_to_save": True, "n_circuit": n_circuit,
                                                         "save_root": save_root, "delta_t": delta_t})
        circuit_list = smart_evolver.run(construct_first=True)
        print(circuit_list[-1])

    return smart_evolver.save_name


if __name__ == '__main__0':
    evolve()
