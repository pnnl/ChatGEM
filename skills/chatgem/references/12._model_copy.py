# Description:
# This script loads a metabolic model from an SBML file (iAF1260.xml) using COBRApy,
# and demonstrates how to create a deep copy of the model object with Python's `copy.deepcopy`.
# This allows you to safely modify `model_copy` without affecting the original model instance.

# Sample ChatGPT prompts to get this code:
# - Python code to make a deep copy of a COBRApy metabolic model.
# - COBRApy: how to copy a metabolic model in memory?
# - How can I duplicate a model in COBRApy for safe editing?
# - Show me COBRApy code to create a copy of a loaded SBML model.
# - Make a deep copy of a COBRApy model using Python.

import os
import cobra
import copy
_ASSETS = os.path.join(os.path.dirname(__file__), '..', 'assets')
def main():
    # Load the original model
    model = cobra.io.read_sbml_model(os.path.join(_ASSETS, 'iAF1260.xml'))
    # Make a deep copy of the model
    model_copy = copy.deepcopy(model)
    print("Model copied. You can now modify model_copy without changing the original model.")

if __name__ == "__main__":
    main()