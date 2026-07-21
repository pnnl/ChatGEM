# Description:
# This script loads the iAF1260 metabolic model from SBML using COBRApy, and uses
# `single_reaction_deletion` to simulate the knockout of a single reaction (PGI).
# The output includes the change in growth rate and essentiality status, useful for
# targeted reaction essentiality analysis in genome-scale models.

# Sample ChatGPT prompts to get this code:
# - Python code to simulate single reaction knockout in COBRApy and print results.
# - COBRApy: run single_reaction_deletion for specific reaction IDs.
# - Example for reaction essentiality analysis in a metabolic model using COBRApy.
# - Print growth rate after single reaction knockout in COBRApy.
# - Script for targeted reaction deletion in COBRApy and show results.

import cobra
from cobra.flux_analysis import single_reaction_deletion
# Load the model
model = cobra.io.read_sbml_model("iAF1260.xml")
# Run single reaction knockout simulation for PGI
deletion_result = single_reaction_deletion(model, reaction_list=["PGI"], processes=1)
print(deletion_result)