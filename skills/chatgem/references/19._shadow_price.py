# Description:
# This script loads a metabolic model from an SBML file (iAF1260.xml) using COBRApy,
# optimizes the model, and saves the shadow prices (dual values for metabolite constraints)
# to a text file called "shadow_prices.txt". Shadow prices help understand how changes in
# metabolite levels affect the objective function.

# Sample ChatGPT prompts to get this code:
# - Python code to optimize a COBRApy model and export shadow prices.
# - COBRApy: how do I save shadow prices to a text file?
# - Show me code to get and export shadow prices from an SBML model using COBRApy.
# - Save the shadow prices from model optimization using COBRApy.
# - COBRApy script to run FBA and save metabolite shadow prices.

import os
import cobra
_ASSETS = os.path.join(os.path.dirname(__file__), '..', 'assets')
def main():
    # Load the model
    model = cobra.io.read_sbml_model(os.path.join(_ASSETS, 'iAF1260.xml'))
    # Optimize the model
    solution = model.optimize()
    # Save shadow prices to a text file
    with open("shadow_prices.txt", "w") as f:
        for met_id, shadow_price in solution.shadow_prices.items():
            f.write(f"{met_id} : {shadow_price}\n")
    print("Shadow prices have been saved to shadow_prices.txt.")

if __name__ == "__main__":
    main()