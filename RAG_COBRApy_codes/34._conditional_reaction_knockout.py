# Description:
# This script loads the iAF1260 E. coli metabolic model in COBRApy, then simulates the knockout
# of a specific reaction (PGI, phosphoglucose isomerase) by setting its flux bounds to zero.
# The knockout is performed temporarily in a context manager and the model is re-optimized to print the new growth rate.
# This demonstrates reaction essentiality analysis in genome-scale models.

# Sample ChatGPT prompts to get this code:
# - Python code to simulate single reaction knockout and print changed growth rate in COBRApy.
# - COBRApy: how can I knock out a reaction and rerun FBA?
# - Example of conditional reaction knockout and FBA analysis in a metabolic model.
# - COBRApy context manager for reaction deletion and growth effect.
# - Print growth rate after reaction knockout for a COBRApy SBML model.

import cobra
# Load model
model = cobra.io.read_sbml_model("iAF1260.xml")
# Pick a reaction to knock out (PGI = phosphoglucose isomerase in glycolysis)
rxn_id = "PGI"
rxn = model.reactions.get_by_id(rxn_id)
# Perform knockout in a context manager so it's temporary
with model:
    rxn.knock_out()
    sol = model.optimize()
    print(f"Growth after knockout of {rxn_id}: {sol.objective_value}")