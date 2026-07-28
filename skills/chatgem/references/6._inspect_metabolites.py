# Description:
# This script loads an SBML metabolic model (iAF1260.xml) using COBRApy and retrieves details
# of a specified metabolite (by its ID). It prints the metabolite's ID, name, formula, and all
# reactions in which it participates. This is useful for investigating metabolite roles and connectivity
# in genome-scale metabolic models.

# Sample ChatGPT prompts to get this code:
# - Python code to get metabolite details and associated reactions in COBRApy.
# - COBRApy: how to print metabolite info and list its reactions?
# - Show COBRApy code to list a metabolite's properties and reactions.
# - How can I display reactions involving a metabolite in COBRApy?
# - Get formula, name, and reactions of a metabolite from SBML using COBRApy.

import os
import cobra
_ASSETS = os.path.join(os.path.dirname(__file__), '..', 'assets')
def main():
    # Load the model
    model = cobra.io.read_sbml_model(os.path.join(_ASSETS, 'iAF1260.xml'))
    # Get the metabolite by its ID
    met_id = "glc__D_c"  # Change this if you want a different metabolite
    met = model.metabolites.get_by_id(met_id)
    # Print metabolite information
    print("ID:", met.id)
    print("Name:", met.name)
    print("Formula:", met.formula)
    print("Reactions:", [r.id for r in met.reactions])

if __name__ == "__main__":
    main()