# Description:
# This script loads a metabolic model from an SBML file ("iAF1260.xml") using COBRApy,
# runs standard Flux Balance Analysis (FBA) and parsimonious FBA (pFBA), then writes a comparison
# of the resulting flux distributions for every reaction to "flux_comparison.txt".
# This allows assessment of differences between typical and parsimonious flux solutions.

# Sample ChatGPT prompts to get this code:
# - Python code to compare FBA and pFBA fluxes in a COBRApy model and save results.
# - COBRApy: how do I run FBA and pFBA and write fluxes to a file?
# - Compare standard and parsimonious FBA flux distributions using COBRApy.
# - Write FBA vs pFBA flux results for all reactions with COBRApy.
# - COBRApy example to save FBA/pFBA reaction fluxes from SBML to text.

import cobra
from cobra.flux_analysis import pfba
def main():
    # Load the model
    model = cobra.io.read_sbml_model("iAF1260.xml")
    # Run FBA
    fba_sol = model.optimize()
    # Run pFBA
    pfba_sol = pfba(model)
    # Write flux comparison to file
    with open("flux_comparison.txt", "w") as f:
        f.write("Reaction\tFBA_flux\tpFBA_flux\n")
        for r in model.reactions:
            fba_flux = fba_sol.fluxes[r.id] if r.id in fba_sol.fluxes else "NA"
            pfba_flux = pfba_sol.fluxes[r.id] if r.id in pfba_sol.fluxes else "NA"
            f.write(f"{r.id}\t{fba_flux}\t{pfba_flux}\n")
    print("FBA and pFBA flux distributions compared and saved to flux_comparison.txt.")

if __name__ == "__main__":
    main()