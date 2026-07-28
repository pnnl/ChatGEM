# Description:
# This script loads a metabolic model from an SBML file (iAF1260.xml) using COBRApy.
# It writes each reaction's ID and its associated subsystem (pathway or functional group)
# to a UTF-8 encoded text file "reaction_subsystems.txt".
# This helps you organize and analyze reactions according to their metabolic subsystems.

# Sample ChatGPT prompts to get this code:
# - Python code to export reactions and their subsystems from a COBRApy model.
# - COBRApy: list all reactions and their subsystems to a text file.
# - How can I write reaction and subsystem information from SBML using COBRApy?
# - Save reaction IDs and subsystems with COBRApy to a file.
# - Export metabolic model reaction subsystems with Python and COBRApy.

import os
import cobra
_ASSETS = os.path.join(os.path.dirname(__file__), '..', 'assets')
def main():
    # Load the model
    model = cobra.io.read_sbml_model(os.path.join(_ASSETS, 'iAF1260.xml'))
    # Write reactions and their subsystems to a text file (UTF-8 encoding)
    with open("reaction_subsystems.txt", "w", encoding="utf-8") as f:
        for r in model.reactions:
            f.write(f"{r.id} → {r.subsystem}\n")
    print("Reaction subsystems have been saved to reaction_subsystems.txt.")

if __name__ == "__main__":
    main()