# Description:
# This script loads the iAF1260 E. coli metabolic model with COBRApy, checks for and adds missing metabolites
# (glyceraldehyde 3-phosphate, G3P) in both extracellular and cytosolic compartments, then adds corresponding 
# exchange and transport reactions (if absent). The updated model is optimized and fluxes for these reactions 
# are printed. This workflow demonstrates how to programmatically expand models and analyze effects of new metabolites/reactions.

# Sample ChatGPT prompts to get this code:
# - Python code to add metabolites and transport/exchange reactions to a COBRApy model, then run FBA.
# - COBRApy: how do I insert new extracellular/cytosolic metabolites and associated reactions?
# - Add metabolite and reactions to a model, optimize, and print new fluxes in COBRApy.
# - Expand a metabolic model with specific transport/exchange steps and show flux distributions.
# - Example: add G3P transport/exchange reactions to SBML model and analyze growth/fba in COBRApy.

import os
import cobra
from cobra import Reaction, Metabolite
_ASSETS = os.path.join(os.path.dirname(__file__), '..', 'assets')
# Load model
model = cobra.io.read_sbml_model(os.path.join(_ASSETS, 'iAF1260.xml'))
# Step 1. Add extracellular metabolite if missing
met_e_id = "g3p_e"
if met_e_id not in model.metabolites:
    met_e = Metabolite(
        met_e_id,
        formula="C3H7O6P",  # Check for proper formula if you wish!
        name="Glyceraldehyde 3-phosphate (extracellular)",
        compartment="e"
    )
    model.add_metabolites([met_e])
    print(f"Added metabolite {met_e_id}")
else:
    met_e = model.metabolites.get_by_id(met_e_id)
    print(f"Metabolite {met_e_id} already present.")
# Step 2. Add exchange reaction (extracellular) if missing
ex_rxn_id = "EX_g3p_e"
if ex_rxn_id not in model.reactions:
    ex_rxn = Reaction(ex_rxn_id)
    ex_rxn.name = "Exchange reaction for G3P (extracellular)"
    ex_rxn.add_metabolites({met_e: -1.0})   # uptake = negative
    ex_rxn.lower_bound = -10.0              # allow import
    ex_rxn.upper_bound = 1000.0             # allow secretion
    model.add_reactions([ex_rxn])
    print(f"Added exchange reaction {ex_rxn_id}")
else:
    print(f"Exchange reaction {ex_rxn_id} already present.")
# Step 3. Add cytosolic metabolite if missing
met_c_id = "g3p_c"
if met_c_id not in model.metabolites:
    met_c = Metabolite(
        met_c_id,
        formula="C3H7O6P",  # Usually same formula as extracellular
        name="Glyceraldehyde 3-phosphate (cytosol)",
        compartment="c"
    )
    model.add_metabolites([met_c])
    print(f"Added metabolite {met_c_id}")
else:
    met_c = model.metabolites.get_by_id(met_c_id)
    print(f"Metabolite {met_c_id} already present.")
# Step 4. Add transport reaction (e → c) if missing
trans_rxn_id = "G3Ptex"
if trans_rxn_id not in model.reactions:
    trans_rxn = Reaction(trans_rxn_id)
    trans_rxn.name = "G3P transport from extracellular to cytosol"
    trans_rxn.add_metabolites({met_e: -1.0, met_c: 1.0})  # e → c
    trans_rxn.lower_bound = 0.0    # Usually only import allowed
    trans_rxn.upper_bound = 1000.0
    model.add_reactions([trans_rxn])
    print(f"Added transport reaction {trans_rxn_id}")
else:
    print(f"Transport reaction {trans_rxn_id} already present.")
# Step 5. Optimize model and show fluxes
solution = model.optimize()
print(f"\nModel status: {solution.status}")
print(f"Objective value (growth rate): {solution.objective_value:.6f}")
# Show fluxes
print(f"Flux through {ex_rxn_id}: {solution.fluxes.get(ex_rxn_id, 0):.6f}")
print(f"Flux through {trans_rxn_id}: {solution.fluxes.get(trans_rxn_id, 0):.6f}")