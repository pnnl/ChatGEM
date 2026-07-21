# Description:
# This script loads the iAF1260 E. coli metabolic model from SBML using COBRApy, and uses
# `double_reaction_deletion` to simulate the double knockout of specified reactions ("PGI" and "PFK").
# The analysis shows the impact of knocking out both reactions (such as growth rate changes),
# useful for studying synthetic lethality and pathway robustness in metabolic networks.

# Sample ChatGPT prompts to get this code:
# - Python code to run double reaction knockout for specified reaction pairs and print results in COBRApy.
# - COBRApy: simulate double reaction deletion and display effects.
# - Example of targeted double reaction knockout analysis in a metabolic model using COBRApy.
# - Print synthetic lethality results for specific reaction pairs in COBRApy.
# - Run double_reaction_deletion for chosen reactions and show output (COBRApy).

import cobra
from cobra.flux_analysis import double_reaction_deletion
# Load model
model = cobra.io.read_sbml_model("iAF1260.xml")
# Double reaction knockout simulation for PGI and PFK
reaction_list1 = ["PGI"]   # First reaction(s) to knockout
reaction_list2 = ["PFK"]   # Second reaction(s) to knockout
# Run double reaction deletion analysis (set processes=1 for Jupyter/Thonny/Windows)
double_rxn_del = double_reaction_deletion(
    model, reaction_list1=reaction_list1, reaction_list2=reaction_list2, processes=1
)
print(double_rxn_del)