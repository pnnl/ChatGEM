# Description:
# This script loads a metabolic model from an SBML file (iAF1260.xml) using COBRApy,
# and prints the model's objective function (as a linear expression) and its direction
# (maximization or minimization). This is useful for understanding the default optimization
# settings of genome-scale metabolic models.

# Sample ChatGPT prompts to get this code:
# - Python code to print the objective function and direction of a COBRApy model.
# - How to display the objective expression and optimization direction in COBRApy?
# - Show me COBRApy code to print metabolic model objective details.
# - Print the SBML model objective function and direction using COBRApy.
# - COBRApy: get and display model objective expression and direction.

import os
import cobra
_ASSETS = os.path.join(os.path.dirname(__file__), '..', 'assets')
def main():
    # Load the model
    model = cobra.io.read_sbml_model(os.path.join(_ASSETS, 'iAF1260.xml'))
    # Print the objective function (linear expression)
    print("Objective expression:", model.objective.expression)
    # Print the direction ('max' or 'min')
    print("Objective direction:", model.objective_direction)

if __name__ == "__main__":
    main()