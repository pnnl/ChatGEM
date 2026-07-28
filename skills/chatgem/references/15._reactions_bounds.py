# Description:
# This script loads a metabolic model from an SBML file (iAF1260.xml) using COBRApy.
# It iterates through all reactions in the model and writes each reaction's ID and flux bounds
# (lower and upper bound) to a text file called "reaction_bounds.txt".
# This helps you quickly reference the allowed flux ranges for every reaction in the model.

# Sample ChatGPT prompts to get this code:
# - Python code to export reaction flux bounds from a COBRApy model to a file.
# - COBRApy: write all reaction lower and upper bounds to a text file.
# - How can I save reaction bounds from a metabolic model using COBRApy?
# - List all reaction flux ranges and write to a file with COBRApy.
# - COBRApy script to write reaction bounds from SBML model to text file.

import os
import cobra
_ASSETS = os.path.join(os.path.dirname(__file__), '..', 'assets')
def main():
    # Load the model
    model = cobra.io.read_sbml_model(os.path.join(_ASSETS, 'iAF1260.xml'))
    # Write reaction bounds to a text file
    with open("reaction_bounds.txt", "w") as f:
        for r in model.reactions:
            f.write(f"{r.id} : {r.lower_bound} < flux < {r.upper_bound}\n")
    print("Reaction bounds have been saved to reaction_bounds.txt.")

if __name__ == "__main__":
    main()