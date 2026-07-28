# Description:
# This script uses COBRApy to load a metabolic model from an SBML file (iAF1260.xml),
# runs an optimization (such as Flux Balance Analysis), and prints the solution status.
# The solution status lets you know whether the optimization succeeded (e.g., 'optimal')
# or failed (e.g., 'infeasible').

# Sample ChatGPT prompts to get this code:
# - Python code to run model optimization and print solution status using COBRApy.
# - COBRApy: how do I check the status of model.optimize()?
# - Show me a script to optimize a COBRApy model and print if it succeeded.
# - Print solution status after FBA with COBRApy.
# - COBRApy SBML loader and optimizer with status output.

import os
import cobra
_ASSETS = os.path.join(os.path.dirname(__file__), '..', 'assets')
def main():
    # Load the model
    model = cobra.io.read_sbml_model(os.path.join(_ASSETS, 'iAF1260.xml'))
    # Optimize the model
    solution = model.optimize()
    # Print solution status
    print("Solution status:", solution.status)

if __name__ == "__main__":
    main()