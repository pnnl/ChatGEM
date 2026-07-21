# Description:
# This script uses COBRApy to load a metabolic model from an SBML file ("iAF1260.xml").
# It exports all reaction IDs and their equations to a text file called "reactions.txt".
# This is useful for quickly listing and reviewing all reactions present in a genome-scale metabolic model.

# Sample ChatGPT prompts to get this code:
# - Python code to export all reaction IDs and equations from a COBRApy model.
# - COBRApy: how to save all reactions and their formulas to a text file?
# - List all reaction equations in a metabolic model using COBRApy.
# - Show COBRApy code to write reaction info from SBML to a file.
# - Save all reactions from a COBRApy model to a text file.

import cobra
def main():
    # Load the model
    model = cobra.io.read_sbml_model("iAF1260.xml")
    # Write all reaction IDs and equations to a text file
    with open("reactions.txt", "w") as f:
        for r in model.reactions:
            f.write(f"{r.id} : {r.reaction}\n")
    print("All reactions have been saved to reactions.txt.")

if __name__ == "__main__":
    main()