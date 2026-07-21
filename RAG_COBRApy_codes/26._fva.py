# Description:
# This script loads a metabolic model from an SBML file ("iAF1260.xml") using COBRApy,
# performs Flux Variability Analysis (FVA) at 100% of the optimum (fraction_of_optimum=1.0),
# prints the first few rows of the FVA result, and saves the full result to a tab-delimited text file.
# This workflow helps identify minimum and maximum possible fluxes for each reaction under optimal growth.

# Sample ChatGPT prompts to get this code:
# - Python code to run FVA in COBRApy and save results to a file.
# - COBRApy: perform flux variability analysis and export to text.
# - How do I perform and save FVA results for a metabolic model in COBRApy?
# - Example script for running FVA and exporting output from SBML using COBRApy.
# - Save flux variability analysis from a COBRApy model to a .txt file.

import cobra
from cobra.flux_analysis import flux_variability_analysis
def main():
    # Load the model
    model = cobra.io.read_sbml_model("iAF1260.xml")
    # Perform FVA for fraction_of_optimum=1.0 (full optimum)
    fva_results = flux_variability_analysis(model, fraction_of_optimum=1.0)
    print(fva_results.head())
    # Save FVA results to a text file
    fva_results.to_csv("fva_results.txt", sep="\t")
    print("FVA results saved to fva_results.txt.")

if __name__ == "__main__":
    main()