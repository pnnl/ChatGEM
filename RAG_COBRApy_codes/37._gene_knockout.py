# Description:
# This script loads the iAF1260 E. coli metabolic model using COBRApy, simulates the knockout
# of a specified gene (by BiGG gene ID) using the `single_gene_deletion` function, and prints
# out the resulting growth rate and status. This is useful for targeted gene essentiality analysis
# in genome-scale models.

# Sample ChatGPT prompts to get this code:
# - Python code to simulate knockout of a specified gene and print results in COBRApy.
# - COBRApy: run single_gene_deletion for a list of gene IDs and show output.
# - How can I test a gene knockout's effect on growth using COBRApy?
# - Print growth rate after specific gene knockout in COBRApy.
# - Example for targeted gene deletion analysis in COBRApy and print results.

import cobra
from cobra.flux_analysis import single_gene_deletion
# Load model
model = cobra.io.read_sbml_model("iAF1260.xml")
# Simulate knockout of a single gene (by BiGG gene ID)
gene_list = ["b3959"]  # Replace with any gene ID of interest
deletion_result = single_gene_deletion(model, gene_list=gene_list, processes=1)
print(deletion_result)