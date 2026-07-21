# Description:
# This script loads a metabolic model from an SBML file ("iAF1260.xml") using the COBRApy library.
# It writes the IDs and names of all metabolites in the model to a text file called "metabolites.txt".
# This is useful for quickly listing and referring to all metabolites present in a genome-scale metabolic model.

# Sample ChatGPT prompts to get this code:
# - Python code to export all metabolite IDs and names from a COBRApy model.
# - COBRApy: how to get and save all metabolites from a metabolic model?
# - Show me COBRApy code to list and write metabolite info from SBML.
# - Write all metabolites to a file from a COBRApy model.
# - Save metabolite IDs and names from a COBRApy model to a text file.

import cobra
def main():
    # Load the model
    model = cobra.io.read_sbml_model("iAF1260.xml")
    # Write all metabolite IDs and names to a text file
    with open("metabolites.txt", "w") as f:
        for m in model.metabolites:
            f.write(f"{m.id} : {m.name}\n")
    print("All metabolites have been saved to metabolites.txt.")

if __name__ == "__main__":
    main()