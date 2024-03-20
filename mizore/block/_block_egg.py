from importlib import import_module

def get_block_egg(block):
    egg_dict=dict()
    egg_dict["ty"]=block.__class__.__name__
    egg_dict["ar"]=block.get_reconstruct_args()
    egg_dict["pa"]=block.parameter
    egg_dict["in"]=block.is_inverse
    return egg_dict

public_block_module = import_module("mizore.block")

def spawn_block_egg(egg_dict,module_name=None):
    if module_name is None:
        block_module=public_block_module
    else:
        block_module=import_module(module_name)
    block_cls = getattr(block_module,egg_dict["ty"])
    block=block_cls(*egg_dict["ar"])
    block.parameter=egg_dict["pa"]
    block.is_inverse=egg_dict["in"]
    return block

