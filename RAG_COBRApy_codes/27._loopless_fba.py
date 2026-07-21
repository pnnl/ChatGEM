# Description:
# This script loads a metabolic model from an SBML file ("iAF1260.xml") using COBRApy,
# runs standard Flux Balance Analysis (FBA) and loopless FBA, then compares the fluxes
# from both solutions in a pandas DataFrame. It prints the first few lines of the comparison
# and saves the full results to a tab-delimited text file. This allows you to see how enforcing
# loopless flux distributions changes the solution.

# Sample ChatGPT prompts to get this code:
# - Python code to compare FBA and loopless FBA fluxes in COBRApy and save results.
# - COBRApy: how can I run loopless FBA and export the results?
# - Compare and write standard vs loopless FBA solutions for a model with COBRApy.
# - Example script for FBA and loopless FBA comparison using COBRApy and pandas.
# - Save standard and loopless FBA flux distributions to file with COBRApy.

import cobra
import pandas
# Load the model from SBML file
model = cobra.io.read_sbml_model("iAF1260.xml")
# Standard FBA
nominal = model.optimize()
# Loopless FBA (returns a Solution object)
from cobra.flux_analysis.loopless import loopless_solution
loopless = loopless_solution(model)
# Create a DataFrame comparing fluxes
df = pandas.DataFrame({'FBA': nominal.fluxes, 'Loopless FBA': loopless.fluxes})
# Print first few lines
print(df.head())
# Optionally, save to a text file
df.to_csv("fba_vs_loopless_fba.txt", sep="\t")
print("Saved flux comparison to fba_vs_loopless_fba.txt.")