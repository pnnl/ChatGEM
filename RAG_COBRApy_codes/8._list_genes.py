# Description:
# This script uses COBRApy to load a metabolic model from an SBML file ("iAF1260.xml"),
# then extracts all gene IDs from the model and writes them to a text file called "genes.txt".
# This is useful for quickly listing all genes in the model for genome-scale or constraint-based
# analysis.

# Sample ChatGPT prompts to get this code:
# - Python code to export all gene IDs from a COBRApy model to a text file.
# - COBRApy: how to get and save all genes from a metabolic model?
# - Show me COBRApy code to list and write gene IDs from SBML.
# - How can I print all gene IDs in a model using COBRApy?
# - Save all genes from a COBRApy model to a file.

import cobra
def main():
    # Load the model
    model = cobra.io.read_sbml_model("iAF1260.xml")
    # Write all gene IDs to a text file
    with open("genes.txt", "w") as f:
        for g in model.genes:
            f.write(g.id + "\n")
    print("All gene IDs have been saved to genes.txt.")

if __name__ == "__main__":
    main()