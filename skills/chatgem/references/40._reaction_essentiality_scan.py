# Description:
# This script loads the iAF1260 E. coli metabolic model in COBRApy and performs single reaction deletion analysis.
# It compares each reaction knockout's growth rate to wild-type and labels reactions as essential if knockout
# growth is less than 10% of the original. Results (reaction ID, growth rate, essentiality status) are saved
# to "essential_reactions.txt". This helps identify essential reactions for growth in genome-scale models.

# Sample ChatGPT prompts to get this code:
# - Python code to run single reaction knockout essentiality scan and save results in COBRApy.
# - COBRApy: use single_reaction_deletion and write growth/essentiality to a text file.
# - Identify essential reactions by FBA knockout simulation with COBRApy and export results.
# - Essentiality scan: reactions with growth < 10% of WT in COBRApy model.
# - How can I export essential reactions from knockout analysis in COBRApy?

import os
import cobra
from cobra.flux_analysis import single_reaction_deletion
_ASSETS = os.path.join(os.path.dirname(__file__), '..', 'assets')
model = cobra.io.read_sbml_model(os.path.join(_ASSETS, 'iAF1260.xml'))
wt_growth = model.optimize().objective_value
essentiality = single_reaction_deletion(model, processes=1)
threshold = 0.1 * wt_growth
with open("essential_reactions.txt", "w") as f:
    f.write("Reaction\tGrowth rate\tEssentiality\n")
    for _, row in essentiality.iterrows():
        # 'ids' is a set, extract the reaction ID
        rxn_set = row["ids"]
        rxn_id = next(iter(rxn_set)) if isinstance(rxn_set, set) else rxn_set
        essen = "yes" if row['growth'] < threshold else "no"
        f.write(f"{rxn_id}\t{row['growth']:.6f}\t{essen}\n")
print("Results saved to essential_reactions.txt")