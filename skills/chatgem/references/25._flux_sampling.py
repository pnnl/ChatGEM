# Description:
# This script loads a metabolic model from an SBML file (iAF1260.xml) using COBRApy,
# and performs flux sampling using two methods: OptGP and ACHR, each generating 1000 samples.
# It prints the first few samples and saves all samples from each method to separate tab-delimited text files.
# This allows exploration of the solution space of feasible flux distributions in genome-scale models.

# Sample ChatGPT prompts to get this code:
# - Python code to perform ACHR and OptGP flux sampling in a COBRApy model and save results.
# - COBRApy: run flux sampling and save samples to a file.
# - Demonstrate how to sample the solution space with OptGP and ACHR using COBRApy.
# - Example COBRApy script to explore and export sampling results from a model.
# - Save flux sample distributions (OptGP, ACHR) with COBRApy to text files.

import os
import cobra
from cobra.sampling import sample
_ASSETS = os.path.join(os.path.dirname(__file__), '..', 'assets')
def main():
    # Load the model
    model = cobra.io.read_sbml_model(os.path.join(_ASSETS, 'iAF1260.xml'))
    # OptGP sampling
    print("Sampling with OptGP ...")
    samples_optgp = sample(model, n=1000, method="optgp")
    print("First 5 OptGP samples:")
    print(samples_optgp.head())
    # Save to text file
    samples_optgp.to_csv("optgp_samples.txt", sep="\t", index=False)
    print("OptGP samples saved to optgp_samples.txt")
    # ACHR sampling
    print("\nSampling with ACHR ...")
    samples_achr = sample(model, n=1000, method="achr")
    print("First 5 ACHR samples:")
    print(samples_achr.head())
    # Save to text file
    samples_achr.to_csv("achr_samples.txt", sep="\t", index=False)
    print("ACHR samples saved to achr_samples.txt")

if __name__ == "__main__":
    main()