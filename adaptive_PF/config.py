import os

# set the root of the output as the directory of this file
output_root = os.path.dirname(os.path.abspath(__file__))
adaptive_pf_output = output_root + "/adaptive_evolution"
os.environ["adaptive_pf_output"] = adaptive_pf_output