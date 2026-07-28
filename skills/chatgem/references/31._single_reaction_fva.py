# Description:
# This script loads a metabolic model from an SBML file (iAF1260.xml) using COBRApy,
# then performs Flux Variability Analysis (FVA) for a single reaction (PGI).
# The results are printed and also saved to a tab-delimited text file.
# This helps you easily investigate the allowable flux range for a specific reaction.

# Sample ChatGPT prompts to get this code:
# - Python code to do FVA for a single reaction in COBRApy and save results.
# - COBRApy: perform flux variability analysis on one reaction.
# - How do I run FVA for a specific reaction and write output to a file with COBRApy?
# - COBRApy example for single-reaction FVA and export.
# - Output FVA results for one reaction in COBRApy to a text file.

import os
import cobra
from cobra.flux_analysis import flux_variability_analysis
_ASSETS = os.path.join(os.path.dirname(__file__), '..', 'assets')
# Load model
model = cobra.io.read_sbml_model(os.path.join(_ASSETS, 'iAF1260.xml'))
# Perform FVA for a single reaction (PGI)
reaction_id = "PGI"
fva_single = flux_variability_analysis(model, reaction_list=[reaction_id])
print(fva_single)
# Save results to file
fva_single.to_csv("fva_PGI.txt", sep="\t")
print('FVA results for PGI saved to fva_PGI.txt.')