# Description:
# This script loads a metabolic model from an SBML file ("iAF1260.xml") using the COBRApy package.
# After loading, it prints an overview with the model's ID, total number of reactions, metabolites,
# and genes. This is a useful quick summary to confirm the successful load and inspect model content.

# Sample ChatGPT prompts to get this code:
# - Python code to load a COBRApy SBML model and print a summary.
# - COBRApy: show model ID, number of reactions, metabolites, and genes.
# - How can I quickly summarize a metabolic model's contents in COBRApy?
# - Print basic statistics about a loaded COBRApy model.
# - COBRApy script for SBML model summary output.

import cobra
def main():
    # Load the model from SBML file
    model = cobra.io.read_sbml_model("iAF1260.xml")
    print(f"Model loaded: {model.id}, {len(model.reactions)} reactions, {len(model.metabolites)} metabolites, {len(model.genes)} genes")

if __name__ == "__main__":
    main()