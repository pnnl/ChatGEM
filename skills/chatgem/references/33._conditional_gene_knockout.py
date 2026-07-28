# Description:
# This script loads a metabolic model from an SBML file (iAF1260.xml) using COBRApy,
# and simulates the knockout of a specific gene by setting its associated reactions' bounds to zero
# (using a context manager for a temporary knockout). The resulting growth rate after knockout 
# is printed. This workflow demonstrates gene essentiality analysis in genome-scale models.

# Sample ChatGPT prompts to get this code:
# - Python code to simulate single gene knockout and print new growth rate in COBRApy.
# - COBRApy: how do I knock out a gene and rerun FBA?
# - Example of gene deletion and growth analysis in a metabolic model with COBRApy.
# - How to conditionally knock out a gene and optimize in COBRApy?
# - COBRApy code to perform gene essentiality analysis by gene knockout.

import os
import cobra
_ASSETS = os.path.join(os.path.dirname(__file__), '..', 'assets')
# Load model
model = cobra.io.read_sbml_model(os.path.join(_ASSETS, 'iAF1260.xml'))
# Pick a gene to knock out
gene_id = "b1773"  # Replace with any valid gene in your model
# Perform knockout in a context manager so it's temporary
with model:
    gene = model.genes.get_by_id(gene_id)
    gene.knock_out()
    sol = model.optimize()
    print(f"Growth after knocking out {gene_id}: {sol.objective_value}")