import math, numpy
from collections import Iterable


class BlockPool:
    """
    Base class of block pools
    Generate a block block_pool by Args when constructed.
    Attributes:
        blocks: Set of blocks in the block_pool
    Overloads:
        pool1+pool2=pool3 is equivalent to merge the two pools
        block_pool+=block is equivalent to add the block to the block_pool
    """

    def __init__(self, block_iter_or_init_block=None):
        """
        Args:
            block_iter_or_init_block: can be a iterator that yields blocks or a block
        """
        self.blocks = set()
        if isinstance(block_iter_or_init_block, Iterable):
            for block in block_iter_or_init_block:
                self.blocks.add(block)
        else:
            if block_iter_or_init_block != None:
                self.blocks.add(block_iter_or_init_block)
        return

    def __iter__(self):
        return iter(self.blocks)

    def __len__(self):
        return len(self.blocks)

    def __str__(self):
        info = "Pool size:" + str(len(self.blocks)) + "\n"
        for block in self.blocks:
            info += str(block) + "\n"
        return info

    def __iadd__(self, block):
        self.blocks.add(block)
        return self

    def __add__(self, another_pool):
        return BlockPool(set.union(another_pool.blocks, self.blocks))

    def generate_random_reduced_pool(self, n_block=0, percent=0):
        """
        Method to randomly extract blocks in the block block_pool to form a new block_pool
        Args:
            n_block: number of blocks in the new block_pool
            percent: the relative size of the new block_pool with respect to the original block_pool
            Only one of the two arguments needs to be provided. If both provided, n_block will be used with priority 
        """
        if (n_block <= 0 or n_block >= len(self.blocks)) and (percent >= 1 or percent <= 0):
            print("Invalid parameter for random block_pool reduce! Self returned.")
            return self
        if n_block == 0:
            n_block = math.ceil(len(self.blocks) * percent)

        new_block_set = random_select(self.blocks, n_block)
        pool = BlockPool()
        pool.blocks = new_block_set
        return pool

    def merge_with_another_pool(self, another_pool):
        self.blocks = set.union(another_pool.blocks, self.blocks)


def random_select(collection=set(), num=1):
    return set(numpy.random.choice(numpy.array(list(collection), dtype=object), size=num, replace=False))
