# Description:
# This Python script uses the COBRApy library to load a metabolic model from an SBML (.xml) file,
# identify the blocked reactions (those with zero flux in all steady-state solutions),
# and export their IDs to a text file called "blocked_reactions.txt".
# This is useful for analyzing and curating genome-scale metabolic models.

# Sample ChatGPT prompts to get this code:
# - Python script to find blocked reactions in a COBRApy model and save them to a file.
# - COBRApy code to load an SBML model and list blocked reactions in a text file.
# - Can you provide code to identify blocked reactions using COBRApy and export them?
# - Give me a Python example for detecting blocked reactions in a metabolic model.
# - How to use COBRApy to find and save blocked reactions from an SBML model?

import cobra
from cobra.flux_analysis import find_blocked_reactions

def main():
    # Load the model from SBML (.xml) file
    model = cobra.io.read_sbml_model("iAF1260.xml")
    # Find blocked reactions (zero flux in any steady-state solution)
    blocked = find_blocked_reactions(model)
    # Write blocked reactions to a text file
    with open("blocked_reactions.txt", "w") as f:
        for rxn_id in blocked:
            f.write(rxn_id + "\n")
    print("Blocked reactions have been saved to blocked_reactions.txt.")

if __name__ == "__main__":
    main()