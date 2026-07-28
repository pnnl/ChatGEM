# Description:
# This script loads a metabolic model from an SBML file (iAF1260.xml) using COBRApy and retrieves
# detailed information about a specified reaction (by its ID). It prints the reaction's ID, name,
# equation, bounds, subsystem, and associated genes. This helps users quickly inspect reaction
# properties and their genetic context in genome-scale models.

# Sample ChatGPT prompts to get this code:
# - Python code to get reaction details and associated genes in COBRApy.
# - COBRApy: how to print reaction info including bounds and genes?
# - Show me COBRApy code to describe a reaction from a metabolic model.
# - How can I show reaction and gene info for a reaction in COBRApy?
# - Get reaction equation, subsystem, and gene associations from SBML with COBRApy.

import os
import cobra
_ASSETS = os.path.join(os.path.dirname(__file__), '..', 'assets')
def main():
    # Load the model
    model = cobra.io.read_sbml_model(os.path.join(_ASSETS, 'iAF1260.xml'))
    # Get the reaction by its ID
    rxn_id = "PGI"  # Change to any reaction ID you want to inspect
    rxn = model.reactions.get_by_id(rxn_id)
    # Print reaction information
    print(rxn.id, rxn.name)
    print("Equation:", rxn.reaction)
    print("Lower bound:", rxn.lower_bound, "Upper bound:", rxn.upper_bound)
    print("Subsystem:", rxn.subsystem)
    print("Associated genes:", [g.id for g in rxn.genes])

if __name__ == "__main__":
    main()