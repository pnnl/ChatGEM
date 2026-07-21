# Description:
# This script loads the iAF1260 E. coli metabolic model in COBRApy, and uses the 
# `double_gene_deletion` function to simulate knockout of specified gene pairs 
# (here: "b3959" and "b3844"). The result, which includes the growth rate after both genes are 
# knocked out, is printed. This is useful for investigating synthetic lethality and double knockout effects.

# Sample ChatGPT prompts to get this code:
# - Python code to run double gene knockout for specified gene pairs and print results in COBRApy.
# - COBRApy: how do I simulate double gene deletion and show outcomes?
# - Example of targeted double gene knockout analysis in a metabolic model using COBRApy.
# - Print synthetic lethality results for specific gene pairs in COBRApy.
# - Run double_gene_deletion for chosen genes and show output (COBRApy).

import cobra
from cobra.flux_analysis import double_gene_deletion
# Load model
model = cobra.io.read_sbml_model("iAF1260.xml")
# Double gene knockout simulation
gene_list1 = ["b3959"]  # First gene(s) to knockout
gene_list2 = ["b3844"]  # Second gene(s) to knockout
# Run the double gene deletion analysis (set processes=1 for Windows/Jupyter)
double_del = double_gene_deletion(model, gene_list1=gene_list1, gene_list2=gene_list2, processes=1)
print(double_del)