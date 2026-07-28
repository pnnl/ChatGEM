# Description:
# This script loads a metabolic model from an SBML file (iAF1260.xml) using COBRApy.
# It identifies transport reactions, defined as reactions involving metabolites from more than
# one compartment. The script saves all such transport reactions (their IDs and equations) to
# "transport_reactions.txt". This helps analyze movement of compounds between cellular compartments.

# Sample ChatGPT prompts to get this code:
# - Python code to extract transport reactions from a COBRApy model and save them.
# - COBRApy: find and export reactions that involve multiple compartments.
# - How can I list and write transport reactions from SBML model using COBRApy?
# - Identify and save multi-compartment reactions in a COBRApy model.
# - Export transport reaction IDs and equations from a metabolic model in Python.

import os
import cobra
_ASSETS = os.path.join(os.path.dirname(__file__), '..', 'assets')
def main():
    # Load the model
    model = cobra.io.read_sbml_model(os.path.join(_ASSETS, 'iAF1260.xml'))
    # Identify transport reactions (involving >1 compartment)
    transport_rxns = [r for r in model.reactions if len(r.compartments) > 1]
    # Write transport reactions to a text file
    with open("transport_reactions.txt", "w") as f:
        for r in transport_rxns:
            f.write(f"{r.id} : {r.reaction}\n")
    print(f"Transport reactions have been saved to transport_reactions.txt ({len(transport_rxns)} found).")

if __name__ == "__main__":
    main()