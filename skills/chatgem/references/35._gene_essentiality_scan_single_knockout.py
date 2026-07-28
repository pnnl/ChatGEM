# Description:
# This script loads the iAF1260 E. coli model in COBRApy, performs single gene knockout analysis using
# the `single_gene_deletion` function, and compares each knockout strain's growth rate to the wild-type.
# Genes with <10% of wild-type growth are considered essential. The results are saved in 
# "essential_genes.txt" with gene ID, growth rate, and essentiality status. This helps identify
# essential genes in genome-scale metabolic models.

# Sample ChatGPT prompts to get this code:
# - Python code to perform single gene knockout analysis and define essential genes in COBRApy.
# - COBRApy: run single_gene_deletion, compare growth, and write results to file.
# - Identify essential genes by knockout simulation with COBRApy and save results.
# - Script to run gene essentiality analysis and export data using COBRApy.
# - Essentiality scan: genes with growth less than 10% of WT in COBRApy model.

import os
import cobra
from cobra.flux_analysis import single_gene_deletion
_ASSETS = os.path.join(os.path.dirname(__file__), '..', 'assets')
def main():
    # Load model
    model = cobra.io.read_sbml_model(os.path.join(_ASSETS, 'iAF1260.xml'))
    # Wild-type growth rate
    wt_sol = model.optimize()
    wt_growth = wt_sol.objective_value
    # Run single gene knockout scan (set processes=1 if using IDE/Jupyter)
    essentiality = single_gene_deletion(model, processes=1)
    # Determine essentiality (growth <10% of WT)
    threshold = 0.1 * wt_growth
    essentiality['essentiality'] = [
        "yes" if gr < threshold else "no" for gr in essentiality["growth"]
    ]
    # Write to file with 3 columns
    with open("essential_genes.txt", "w") as f:
        f.write("gene\tgrowth_rate\tessentiality\n")
        for gene, row in essentiality.iterrows():
            f.write(f"{gene}\t{row['growth']:.6f}\t{row['essentiality']}\n")
    print("Results saved to essential_genes.txt")

if __name__ == "__main__":
    main()