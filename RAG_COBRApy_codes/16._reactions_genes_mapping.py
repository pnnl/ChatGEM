# Description:
# This script loads a metabolic model from an SBML file ("iAF1260.xml") using COBRApy.
# It writes a mapping of all reactions to their associated gene IDs to a text file 
# named "reaction_gene_mapping.txt". This mapping is useful for analyzing the genetic 
# basis of each reaction within a genome-scale metabolic model.

# Sample ChatGPT prompts to get this code:
# - Python code to export reaction-gene associations from a COBRApy model.
# - COBRApy: write all reaction to gene mappings to a text file.
# - List reaction and their genes for a metabolic model using COBRApy.
# - Save reactions and associated gene IDs from SBML to a text file with COBRApy.
# - COBRApy script to output reaction-gene mapping from a model.

import cobra
def main():
    # Load the model
    model = cobra.io.read_sbml_model("iAF1260.xml")
    # Write reaction-gene mappings to a text file
    with open("reaction_gene_mapping.txt", "w") as f:
        for r in model.reactions:
            gene_ids = [g.id for g in r.genes]
            f.write(f"{r.id} : {gene_ids}\n")
    print("Reaction-Gene mapping has been saved to reaction_gene_mapping.txt.")

if __name__ == "__main__":
    main()