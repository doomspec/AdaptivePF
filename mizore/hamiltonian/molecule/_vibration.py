from openfermion.ops import FermionOperator, QubitOperator
from openfermion.transforms import bravyi_kitaev, jordan_wigner
from openfermion.utils import hermitian_conjugated

from pyscf import gto, dft
from pyscf.prop.freq import rks


def mol_vibration():
    mol = gto.M(atom='''
                O 0 0      0
                H 0 -0.757 0.587
                H 0  0.757 0.587''',
                basis='ccpvdz', verbose=4)
    mf = dft.RKS(mol).run()
    # w:9 list, nodes:3*3*9 ndarray
    w, modes = rks.Freq(mf).kernel()
    return w, modes


def second_order_vibration():  # TODO
    """ 
    many body expansion the vibrational potential, 
    and truncate it to second order for practice 
    """
    hamiltonian = FermionOperator('0^ 3', 0.5)

    return hamiltonian
