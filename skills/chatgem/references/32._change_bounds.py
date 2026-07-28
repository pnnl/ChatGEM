# Description:
# This script loads the iAF1260 E. coli metabolic model using COBRApy, sets new bounds for the glucose exchange reaction
# (limiting import to 5 mmol/gDW/hr and preventing secretion), then runs flux balance analysis (FBA) and prints the resulting growth rate.
# This workflow demonstrates how changing substrate uptake constraints impacts simulated growth in genome-scale models.

# Sample ChatGPT prompts to get this code:
# - Python code to change glucose uptake bounds and optimize a COBRApy model.
# - COBRApy: how do I set exchange reaction limits and rerun FBA?
# - Example of constraining substrate import and showing effect on growth using COBRApy.
# - Change reaction bounds, optimize, and print objective value for SBML model in COBRApy.
# - COBRApy: set exchange lower/upper bound and print new growth rate after FBA.

import os
import cobra
_ASSETS = os.path.join(os.path.dirname(__file__), '..', 'assets')
# Load model
model = cobra.io.read_sbml_model(os.path.join(_ASSETS, 'iAF1260.xml'))
# Change bounds for glucose exchange reaction
rxn = model.reactions.get_by_id("EX_glc__D_e")
rxn.lower_bound = -5   # Limit glucose import to 5 mmol/gDW/hr
rxn.upper_bound = 0    # No secretion, only import
# Optimize the model (FBA)
solution = model.optimize()
print("Growth after bound change:", solution.objective_value)