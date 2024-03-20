from openfermion import QubitOperator
from mizore.objective import EnergyObjective

def full_connected_transverse_field_ising(n_qubit,field=1,spin_coupling=1):
    spin_coupling_hamiltonian=QubitOperator()
    for i in range(n_qubit):
        for j in range(i+1,n_qubit):
            spin_coupling_hamiltonian+=QubitOperator("Z"+str(i))*QubitOperator("Z"+str(j))
    field_hamitonian=QubitOperator()
    for i in range(n_qubit):
        field_hamitonian+=QubitOperator("X"+str(i))
    hamiltonian=spin_coupling*spin_coupling_hamiltonian+field*field_hamitonian
    energy_obj=EnergyObjective(hamiltonian,n_qubit,init_block=None,obj_info={"system_name":"Ising"})
    return energy_obj