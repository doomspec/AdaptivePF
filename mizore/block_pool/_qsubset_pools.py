from ..utilities.Iterators import iter_qsubset_pauli_of_operator
from ._block_pool import BlockPool
import math


def qsubset2number(qsubset):
    value = 0
    for i in qsubset:
        value += 1 << i
    return value


def number2qsubset(number):
    highest = int(math.log2(number)) + 1
    qsubset = []
    for i in range(0, highest):
        if int((number % (1 << (i + 1))) / (1 << i)) == 1:
            qsubset.append(i)
    return qsubset


def get_operator_qsubset_pool(operator):
    qsubsets = set()
    for qsubset, _pauli in iter_qsubset_pauli_of_operator(operator):
        qsubsets.add(qsubset2number(qsubset))
    return qsubsets


def get_qsubset_pool_reduced_block_pool(block_pool, qsubset_pool):
    subset_number_pool = [qsubset2number(subset) for subset in qsubset_pool]
    reduced_pool = BlockPool()
    for block in block_pool:
        if qsubset2number(block.qsubset) in subset_number_pool:
            reduced_pool += block
    return reduced_pool


def iter_entangler_by_qsubsets(qsubsets, constructor):
    for qsubset in qsubsets:
        if len(qsubset) >= 2:
            yield constructor(qsubset)
