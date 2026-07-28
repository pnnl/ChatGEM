# Description:
# This script uses the COBRApy package to load a metabolic model from an SBML file (iAF1260.xml).
# It extracts all exchange reactions from the model and writes their IDs and reaction formulas
# to a text file called "exchange_reactions.txt". This helps users identify uptake/secretion reactions
# in metabolic networks for further analysis.

# Sample ChatGPT prompts to get this code:
# - Python script to export exchange reactions from COBRApy model to a text file.
# - COBRApy code to list all exchange reactions and save them.
# - Show me how to extract exchange reactions from a metabolic model using COBRApy.
# - Can you provide Python code to get exchange reactions and write them to a file?
# - How to use COBRApy to find and export exchange reactions from an SBML file?

import os
import cobra
_ASSETS = os.path.join(os.path.dirname(__file__), '..', 'assets')
def main():
    # Load the model
    model = cobra.io.read_sbml_model(os.path.join(_ASSETS, 'iAF1260.xml'))
    # Get exchange reactions
    exchange_rxns = [r for r in model.exchanges]
    # Write the exchange reactions to a text file
    with open("exchange_reactions.txt", "w") as f:
        for r in exchange_rxns:
            f.write(f"{r.id} : {r.reaction}\n")
    print("Exchange reactions have been saved to exchange_reactions.txt.")

if __name__ == "__main__":
    main()