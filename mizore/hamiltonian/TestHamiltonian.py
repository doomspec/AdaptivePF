from openfermion import MolecularData
from openfermion.transforms import bravyi_kitaev, get_fermion_operator
from .molecule._mizore_run_pyscf import run_pyscf
from .molecule._geometry_generator import geometry_generator_dict, equilibrium_geometry_dict
from .molecule._generate_HF_operation import get_HF_operator
from ..objective._energy_obj import EnergyObjective
from ..block import HartreeFockInitBlock
from ..utilities.Tools import get_operator_qsubset
from mizore.hamiltonian.transforms.FermionTransform import make_transform_spin_separating, get_parity_transform

NOT_DEFINED = 999999
CHEMICAL_ACCURACY = 0.001

"""
The methods for generating simple molecular and graph theory Hamiltonian for VQE to find ground state energy.

make_example_H2, make_example_LiH, make_example_H2O and make_example_N2 are the main methods
Use default parameter will produce a standard Hamiltonian for benchmarking

make_molecular_energy_obj can be used to generate Hamiltonians in a more expert way
Selecting active space based on *Irrep* is implemented in Mizore based on PySCF
Please refer to the document of PySCF to see how to use the irrep symbols
"""


def make_example_H2(basis="sto-3g",
                    geometry_info=equilibrium_geometry_dict["H2"],
                    fermi_qubit_transform=bravyi_kitaev,
                    is_computed=False):
    return make_molecular_energy_obj(molecule_name="H2", basis=basis, geometry_info=geometry_info,
                                     fermi_qubit_transform=fermi_qubit_transform, is_computed=is_computed)


def make_example_H6(basis="sto-3g",
                    geometry_info=equilibrium_geometry_dict["H6"],
                    fermi_qubit_transform=bravyi_kitaev,
                    is_computed=False):
    return make_molecular_energy_obj(molecule_name="H6", basis=basis, geometry_info=geometry_info,
                                     fermi_qubit_transform=fermi_qubit_transform, is_computed=is_computed)


def make_example_H4(basis="sto-3g",
                    geometry_info=equilibrium_geometry_dict["H4"],
                    fermi_qubit_transform=bravyi_kitaev,
                    is_computed=False):
    return make_molecular_energy_obj(molecule_name="H4", basis=basis, geometry_info=geometry_info,
                                     fermi_qubit_transform=fermi_qubit_transform, is_computed=is_computed)


def make_example_LiH(basis="sto-3g",
                     geometry_info=equilibrium_geometry_dict["LiH"],
                     fermi_qubit_transform=bravyi_kitaev,
                     is_computed=False):
    n_cancel_orbital = 2
    n_frozen_orbital = 1
    cas_irrep_nocc = {'A1': 3}
    cas_irrep_ncore = {'E1x': 0, 'E1y': 0}
    return make_molecular_energy_obj(molecule_name="LiH", basis=basis, geometry_info=geometry_info,
                                     n_cancel_orbital=n_cancel_orbital, n_frozen_orbital=n_frozen_orbital,
                                     cas_irrep_nocc=cas_irrep_nocc, cas_irrep_ncore=cas_irrep_ncore,
                                     fermi_qubit_transform=fermi_qubit_transform, is_computed=is_computed)


def make_example_H2O(basis="6-31g",
                     geometry_info=equilibrium_geometry_dict["H2O"],
                     fermi_qubit_transform=bravyi_kitaev,
                     is_computed=False):
    n_cancel_orbital = 5
    n_frozen_orbital = 3
    cas_irrep_nocc = {'B1': 2, 'A1': 3}
    cas_irrep_ncore = {'B1': 0, 'A1': 2}
    return make_molecular_energy_obj(molecule_name="H2O", basis=basis, geometry_info=geometry_info,
                                     n_cancel_orbital=n_cancel_orbital, n_frozen_orbital=n_frozen_orbital,
                                     cas_irrep_nocc=cas_irrep_nocc, cas_irrep_ncore=cas_irrep_ncore,
                                     fermi_qubit_transform=fermi_qubit_transform, is_computed=is_computed)

def make_example_12_spin_orbital_H2O(basis="6-31g",
                     geometry_info=equilibrium_geometry_dict["H2O"],
                     fermi_qubit_transform=bravyi_kitaev,
                     is_computed=False):
    n_cancel_orbital = 5
    n_frozen_orbital = 2
    cas_irrep_nocc = {'A1': 3, 'B1': 1,  "B2": 2}
    cas_irrep_ncore = {'B1': 0, 'A1': 2}
    return make_molecular_energy_obj(molecule_name="H2O", basis=basis, geometry_info=geometry_info,
                                     n_cancel_orbital=n_cancel_orbital, n_frozen_orbital=n_frozen_orbital,
                                     cas_irrep_nocc=cas_irrep_nocc, cas_irrep_ncore=cas_irrep_ncore,
                                     fermi_qubit_transform=fermi_qubit_transform, is_computed=is_computed)


def make_example_N2(basis="cc-pvdz", geometry_info=equilibrium_geometry_dict["N2"], fermi_qubit_transform=bravyi_kitaev,
                    is_computed=False):
    n_cancel_orbital = 18
    n_frozen_orbital = 2
    return make_molecular_energy_obj(molecule_name="N2", basis=basis, geometry_info=geometry_info,
                                     n_cancel_orbital=n_cancel_orbital, n_frozen_orbital=n_frozen_orbital,
                                     fermi_qubit_transform=fermi_qubit_transform, is_computed=is_computed)


def make_2_qubit_H2(geometry_info=equilibrium_geometry_dict["H2"], is_computed=False):
    from . import get_reduced_energy_obj_with_HF_init
    transform = make_transform_spin_separating(get_parity_transform(4), 4)
    energy_obj = make_example_H2(fermi_qubit_transform=transform, geometry_info=geometry_info, is_computed=is_computed)
    energy_obj = get_reduced_energy_obj_with_HF_init(energy_obj, [1, 3])
    return energy_obj


def make_6_qubit_H2(geometry_info=equilibrium_geometry_dict["H2"], is_computed=False):
    from . import get_reduced_energy_obj_with_HF_init
    transform = make_transform_spin_separating(get_parity_transform(8), 8)
    energy_obj = make_example_H2(basis="6-31g", fermi_qubit_transform=transform, geometry_info=geometry_info,
                                 is_computed=is_computed)
    energy_obj = get_reduced_energy_obj_with_HF_init(energy_obj, [3, 7])
    return energy_obj


def make_molecular_energy_obj(molecule_name, basis="sto-3g", geometry_info=None, n_cancel_orbital=0, n_frozen_orbital=0,
                              cas_irrep_nocc=None, cas_irrep_ncore=None, fermi_qubit_transform=bravyi_kitaev,
                              is_computed=False):
    if geometry_info == None:
        geometry_info = equilibrium_geometry_dict[molecule_name]

    # Get geometry
    if molecule_name not in geometry_generator_dict.keys():
        print("No such example molecule, using default H2 hamiltonian.")
        molecule_name = "H2"

    geometry = geometry_generator_dict[molecule_name](geometry_info)

    # Get fermion Hamiltonian

    multiplicity = 1
    charge = 0
    molecule = MolecularData(geometry, basis, multiplicity, charge, str(
        geometry_info))
    molecule.symmetry = True
    if not is_computed:
        molecule = run_pyscf(molecule, run_fci=1, n_frozen_orbital=n_frozen_orbital,
                             n_cancel_orbital=n_cancel_orbital, cas_irrep_nocc=cas_irrep_nocc,
                             cas_irrep_ncore=cas_irrep_ncore, verbose=False)
    molecule.load()

    active_space_start = n_frozen_orbital
    active_space_stop = molecule.n_orbitals - n_cancel_orbital
    n_active_orb = active_space_stop - active_space_start
    molecule.n_orbitals = n_active_orb
    molecule.n_qubits = n_active_orb * 2
    molecule.n_electrons = molecule.n_electrons - active_space_start * 2

    fermion_hamiltonian = get_fermion_operator(
        molecule.get_molecular_hamiltonian(occupied_indices=molecule.frozen_orbitals,
                                           active_indices=molecule.active_orbitals))

    # Map ferimon Hamiltonian to qubit Hamiltonian
    qubit_hamiltonian = fermi_qubit_transform(fermion_hamiltonian)

    # qubit_electron_operator=fermi_qubit_transform(get_electron_fermion_operator(molecule.n_electrons))
    qubit_electron_operator = get_HF_operator(
        molecule.n_electrons, fermi_qubit_transform)
    # qubit_hamiltonian=get_dressed_operator(qubit_electron_operator,qubit_hamiltonian)

    # Ignore terms in Hamiltonian that close to zero
    qubit_hamiltonian.compress()

    # Set the terminate_energy to be achieving the chemical accuracy
    terminate_energy = molecule.fci_energy + CHEMICAL_ACCURACY
    obj_info = {"start_cost": molecule.hf_energy,
                "terminate_cost": terminate_energy,"system_name":molecule_name}

    init_operator = HartreeFockInitBlock(
        get_operator_qsubset(qubit_electron_operator))

    return EnergyObjective(qubit_hamiltonian, molecule.n_qubits, init_operator, obj_info)
