# Description:
# This script loads a metabolic model from an SBML file ("iAF1260.xml") using COBRApy.
# It extracts all compartments in the model and writes their IDs and names to a text file called "compartments.txt".
# This is useful for summarizing the cellular or subcellular compartments included in a genome-scale metabolic model.

# Sample ChatGPT prompts to get this code:
# - Python code to export all compartment IDs and names from a COBRApy model.
# - COBRApy: how to get and save compartments from a metabolic model?
# - List all compartments in a COBRApy model and write them to a file.
# - Show COBRApy code to write compartment info from SBML to a text file.
# - Save compartment IDs and names from COBRApy SBML model to a file.

import cobra
def main():
    # Load the model
    model = cobra.io.read_sbml_model("iAF1260.xml")
    # Write compartments to a text file
    with open("compartments.txt", "w") as f:
        for cid, cname in model.compartments.items():
            f.write(f"{cid} : {cname}\n")
    print("Compartments have been saved to compartments.txt.")

if __name__ == "__main__":
    main()