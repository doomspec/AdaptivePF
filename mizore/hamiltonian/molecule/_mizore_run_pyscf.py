#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

"""Driver to initialize molecular object from pyscf program."""

# Modified by Zijian Zhang
# This attribute can be saved and loaded


from __future__ import absolute_import

from functools import reduce

import numpy
import pyscf
from pyscf import gto, scf, ao2mo, ci, cc, fci, mp
from pyscf import mcscf  # ZZJ CHANGE

from openfermion import MolecularData
from openfermionpyscf import PyscfMolecularData


def prepare_pyscf_molecule(molecule):
    """
    This function creates and saves a pyscf input file.

    Args:
        molecule: An instance of the MolecularData class.

    Returns:
        pyscf_molecule: A pyscf molecule instance.
    """
    pyscf_molecule = gto.Mole()
    pyscf_molecule.atom = molecule.geometry
    pyscf_molecule.basis = molecule.basis
    pyscf_molecule.spin = molecule.multiplicity - 1
    pyscf_molecule.charge = molecule.charge
    pyscf_molecule.symmetry = molecule.symmetry  # ZZJ
    pyscf_molecule.build(symmetry=molecule.symmetry)
    print("Symmetry:", pyscf_molecule.topgroup,
          " is used when build the molecule.")
    return pyscf_molecule


def compute_scf(pyscf_molecule):
    """
    Perform a Hartree-Fock calculation.

    Args:
        pyscf_molecule: A pyscf molecule instance.

    Returns:
        pyscf_scf: A PySCF "SCF" calculation object.
    """
    if pyscf_molecule.spin:
        pyscf_scf = scf.ROHF(pyscf_molecule)
    else:
        pyscf_scf = scf.RHF(pyscf_molecule)
    return pyscf_scf


def compute_integrals(pyscf_molecule, pyscf_scf):
    """
    Compute the 1-electron and 2-electron integrals.

    Args:
        pyscf_molecule: A pyscf molecule instance.
        pyscf_scf: A PySCF "SCF" calculation object.

    Returns:
        one_electron_integrals: An N by N array storing h_{pq}
        two_electron_integrals: An N by N by N by N array storing h_{pqrs}.
    """
    # Get one electrons integrals.
    n_orbitals = pyscf_scf.mo_coeff.shape[1]
    one_electron_compressed = reduce(numpy.dot, (pyscf_scf.mo_coeff.T,
                                                 pyscf_scf.get_hcore(),
                                                 pyscf_scf.mo_coeff))
    one_electron_integrals = one_electron_compressed.reshape(
        n_orbitals, n_orbitals).astype(float)

    # Get two electron integrals in compressed format.
    two_electron_compressed = ao2mo.kernel(pyscf_molecule,
                                           pyscf_scf.mo_coeff)

    two_electron_integrals = ao2mo.restore(
        1,  # no permutation symmetry
        two_electron_compressed, n_orbitals)
    # See PQRS convention in OpenFermion.hamiltonians._molecular_data
    # h[p,q,r,s] = (ps|qr)
    two_electron_integrals = numpy.asarray(
        two_electron_integrals.transpose(0, 2, 3, 1), order='C')

    # Return.
    return one_electron_integrals, two_electron_integrals


def run_pyscf(molecule,
              run_fci=False,
              verbose=False,
              n_frozen_orbital=0, n_cancel_orbital=0, cas_irrep_nocc=None, cas_irrep_ncore=None):  # ZZJ CHANGE
    """
    This function runs a pyscf calculation.

    Args:
        molecule: An instance of the MolecularData or PyscfMolecularData class.
        run_scf: Optional boolean to run SCF calculation.
        run_mp2: Optional boolean to run MP2 calculation.
        run_cisd: Optional boolean to run CISD calculation.
        run_ccsd: Optional boolean to run CCSD calculation.
        run_fci: Optional boolean to FCI calculation.  ## NOTICE !!! FCI is changed to only done in the active space --- ZZJ CHANGE
        verbose: Boolean whether to print calculation results to screen.

    Returns:
        molecule: The updated PyscfMolecularData object. Note the attributes
        of the input molecule are also updated in this function.
    """
    run_scf = True  # ZZJ CHANGE scf must be run
    # Prepare pyscf molecule.
    pyscf_molecule = prepare_pyscf_molecule(molecule)
    molecule.n_orbitals = int(pyscf_molecule.nao_nr())
    molecule.n_qubits = 2 * molecule.n_orbitals
    molecule.nuclear_repulsion = float(pyscf_molecule.energy_nuc())

    # Run SCF.
    pyscf_scf = compute_scf(pyscf_molecule)
    pyscf_scf.verbose = 0
    pyscf_scf.run()
    molecule.hf_energy = float(pyscf_scf.e_tot)
    if verbose:
        print('Hartree-Fock energy for {} ({} electrons) is {}.'.format(
            molecule.name, molecule.n_electrons, molecule.hf_energy))

    # Hold pyscf data in molecule. They are required to compute density
    # matrices and other quantities.
    molecule._pyscf_data = pyscf_data = {}
    pyscf_data['mol'] = pyscf_molecule
    pyscf_data['scf'] = pyscf_scf

    # Populate fields.
    molecule.canonical_orbitals = pyscf_scf.mo_coeff.astype(float)
    molecule.orbital_energies = pyscf_scf.mo_energy.astype(float)

    # Get integrals.
    one_body_integrals, two_body_integrals = compute_integrals(
        pyscf_molecule, pyscf_scf)
    molecule.one_body_integrals = one_body_integrals
    molecule.two_body_integrals = two_body_integrals
    molecule.overlap_integrals = pyscf_scf.get_ovlp()

    # Run FCI.
    if run_fci:
        # pyscf_fci = fci.FCI(pyscf_scf)
        # pyscf_fci.verbose = 0
        # molecule.fci_energy, molecule.fci_wfn = pyscf_fci.kernel()
        # pyscf_data['fci'] = pyscf_fci
        # pyscf_data['fci_wfn'] = molecule.fci_wfn

        nelec = molecule.n_electrons - n_frozen_orbital * 2
        norb = molecule.n_orbitals - n_cancel_orbital - n_frozen_orbital
        pyscf_fci = mcscf.CASCI(pyscf_scf, norb, nelec)
        pyscf_fci.verbose = 0
        if cas_irrep_nocc != None:
            if verbose:
                orbsym = numpy.asarray(scf.hf_symm.get_orbsym(pyscf_scf.mol, pyscf_scf.mo_coeff))
                print("orbsym",orbsym)
            molecule.active_orbitals = mcscf.caslst_by_irrep(
                pyscf_fci, pyscf_scf.mo_coeff, cas_irrep_nocc, cas_irrep_ncore)
            mo = pyscf_fci.sort_mo_by_irrep(cas_irrep_nocc, cas_irrep_ncore)
            frozen_orbitals_0 = [i for i in range(molecule.active_orbitals[norb - 1] + 1)]
            for i in range(len(molecule.active_orbitals)):
                molecule.active_orbitals[i] -= 1
                frozen_orbitals_0[molecule.active_orbitals[i]] = -1
            molecule.frozen_orbitals = []
            for orb in frozen_orbitals_0:
                if orb != -1:
                    molecule.frozen_orbitals.append(orb)
                if n_frozen_orbital == len(molecule.frozen_orbitals):
                    break

            pyscf_data['active_orbitals'] = molecule.active_orbitals
            pyscf_data['frozen_orbitals'] = molecule.frozen_orbitals
        else:
            mo = None
            molecule.active_orbitals = [i for i in range(
                n_frozen_orbital, n_frozen_orbital + norb)]
            molecule.frozen_orbitals = [i for i in range(n_frozen_orbital)]

        fci_result = pyscf_fci.kernel(mo)
        pyscf_fci.analyze()
        molecule.fci_energy = fci_result[0]
        molecule.fci_wfn = fci_result[2]

        if verbose:
            print('FCI energy for {} ({} electrons) is {}.'.format(
                molecule.name, molecule.n_electrons, molecule.fci_energy))

    # Return updated molecule instance.
    pyscf_molecular_data = PyscfMolecularData.__new__(PyscfMolecularData)
    pyscf_molecular_data.__dict__.update(molecule.__dict__)
    pyscf_molecular_data.save()
    return pyscf_molecular_data


def generate_molecular_hamiltonian(
        geometry,
        basis,
        multiplicity,
        charge=0,
        n_active_electrons=None,
        n_active_orbitals=None):
    """Generate a molecular Hamiltonian with the given properties.

    Args:
        geometry: A list of tuples giving the coordinates of each atom.
            An example is [('H', (0, 0, 0)), ('H', (0, 0, 0.7414))].
            Distances in angstrom. Use atomic symbols to
            specify atoms.
        basis: A string giving the basis set. An example is 'cc-pvtz'.
            Only optional if loading from file.
        multiplicity: An integer giving the spin multiplicity.
        charge: An integer giving the charge.
        n_active_electrons: An optional integer specifying the number of
            electrons desired in the active space.
        n_active_orbitals: An optional integer specifying the number of
            spatial orbitals desired in the active space.

    Returns:
        The Hamiltonian as an InteractionOperator.
    """

    # Run electronic structure calculations
    molecule = run_pyscf(
        MolecularData(geometry, basis, multiplicity, charge)
    )

    # Freeze core orbitals and truncate to active space
    if n_active_electrons is None:
        n_core_orbitals = 0
        occupied_indices = None
    else:
        n_core_orbitals = (molecule.n_electrons - n_active_electrons) // 2
        occupied_indices = list(range(n_core_orbitals))

    if n_active_orbitals is None:
        active_indices = None
    else:
        active_indices = list(range(n_core_orbitals,
                                    n_core_orbitals + n_active_orbitals))

    return molecule.get_molecular_hamiltonian(
        occupied_indices=occupied_indices,
        active_indices=active_indices)
