# Description:
# This script loads the iAF1260 E. coli model in COBRApy, performs double gene knockout analysis
# using the `double_gene_deletion` function over all genes (or a selectable subset), and compares
# each knockout strain's growth rate to wild-type. Pairs with <10% of wild-type growth are labeled
# essential. Results (gene1, gene2, growth rate, essentiality) are saved to "double_gene_essentiality.txt".
# This enables discovery of synthetic lethal gene pairs in genome-scale metabolic models.

# Sample ChatGPT prompts to get this code:
# - Python code for double gene knockout essentiality scan and save results in COBRApy.
# - COBRApy: run double_gene_deletion and write growth/essentiality to a text file.
# - Identify essential double gene knockouts using COBRApy and export all pairs.
# - Script for double knockout and synthetic lethality screen using COBRApy.
# - COBRApy: scan all gene pairs and save double gene essentiality results.

import os
import cobra
from cobra.flux_analysis import double_gene_deletion
_ASSETS = os.path.join(os.path.dirname(__file__), '..', 'assets')
def main():
    # Load model
    model = cobra.io.read_sbml_model(os.path.join(_ASSETS, 'iAF1260.xml'))
    # Wild-type growth rate
    wt_sol = model.optimize()
    wt_growth = wt_sol.objective_value
    # Run double gene knockout scan
    # Get list of all gene IDs
    genes_list = [g.id for g in model.genes]
    # Optionally, use a subset for a faster test (e.g. genes_list[:100])
    # result = double_gene_deletion(model, gene_list=genes_list[:100], processes=1)
    result = double_gene_deletion(model, gene_list=genes_list, processes=1)
    # Determine essentiality (growth <10% of WT)
    threshold = 0.1 * wt_growth
    result['essentiality'] = [
        "yes" if gr < threshold else "no" for gr in result["growth"]
    ]
    # Write to file with four columns
    with open("double_gene_essentiality.txt", "w") as f:
        f.write("gene1\tgene2\tgrowth_rate\tessentiality\n")
        for i, row in result.iterrows():
            f.write(f"{row['gene_1']}\t{row['gene_2']}\t{row['growth']:.6f}\t{row['essentiality']}\n")
    print("Results saved to double_gene_essentiality.txt")

if __name__ == "__main__":
    main()