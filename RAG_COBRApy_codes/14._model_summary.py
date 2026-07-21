# Description:
# This script uses COBRApy to load a metabolic model from an SBML file ("iAF1260.xml").
# It prints key statistics including model ID, numbers of reactions, metabolites, and genes,
# and the objective function expression. It then shows a full model summary using COBRApy's summary method,
# providing an overview of network structure and properties.

# Sample ChatGPT prompts to get this code:
# - Python code to print basic stats and show the summary of a COBRApy model.
# - COBRApy: how to display model statistics and summary after loading SBML?
# - Show me a script to print model details and call the summary method in COBRApy.
# - Print the objective, numbers of genes/metabolites/reactions, and run model.summary() with COBRApy.
# - COBRApy SBML loader with statistics and summary output.

import cobra
def main():
    # Load the model
    model = cobra.io.read_sbml_model("iAF1260.xml")
    print("Model ID:", model.id)
    print("Reactions:", len(model.reactions))
    print("Metabolites:", len(model.metabolites))
    print("Genes:", len(model.genes))
    print("Objective:", model.objective.expression)
    print()  # Blank line for readability
    # Show the model summary
    model.summary()

if __name__ == "__main__":
    main()