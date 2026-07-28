# Description:
# This script uses COBRApy to load a metabolic model from an SBML file (iAF1260.xml).
# It retrieves information about a specific gene (by gene ID), printing the gene's ID
# and listing all reactions associated with that gene. This is useful for identifying the
# metabolic role and context of particular genes in genome-scale models.

# Sample ChatGPT prompts to get this code:
# - Python code to get reactions associated with a gene in COBRApy.
# - COBRApy: how to print reactions related to a specific gene?
# - Show me COBRApy code to display gene info and its reactions from a metabolic model.
# - How can I list reactions linked to a gene in a COBRApy model?
# - Get gene details and related reactions from SBML with COBRApy.

import os
import cobra
_ASSETS = os.path.join(os.path.dirname(__file__), '..', 'assets')
def main():
    # Load the model
    model = cobra.io.read_sbml_model(os.path.join(_ASSETS, 'iAF1260.xml'))
    # Get the gene by its ID
    gene_id = "b3599"  # Change this to your gene of interest
    gene = model.genes.get_by_id(gene_id)
    # Display gene information
    print("Gene ID:", gene.id)
    print("Associated reactions:", [r.id for r in gene.reactions])

if __name__ == "__main__":
    main()