# Description:
# This script loads a metabolic model from an SBML file ("iAF1260.xml") using COBRApy,
# automatically finds the biomass reaction by ID, runs parsimonious FBA (pFBA), and prints the model's 
# growth rate (biomass reaction flux). It then saves all reaction fluxes from pFBA to "pfba_fluxes.txt".
# This enables both convenient simulation and documentation of flux results for further analysis.

# Sample ChatGPT prompts to get this code:
# - Python code to run pFBA, print biomass flux, and save all fluxes from a COBRApy model.
# - COBRApy: find biomass reaction, run pFBA, and export results.
# - Example to get and save pFBA fluxes after automatic biomass reaction detection.
# - How to run pFBA in COBRApy and write fluxes plus growth rate.
# - Save pFBA results and print growth rate for metabolic model using COBRApy.

import cobra
from cobra.flux_analysis import pfba
def main():
    model = cobra.io.read_sbml_model("iAF1260.xml")
    # Find biomass reaction ID
    biomass_id = None
    for r in model.reactions:
        if "BIOMASS" in r.id.upper():
            biomass_id = r.id
            print("Biomass reaction found:", biomass_id)
            break
    if not biomass_id:
        raise ValueError("No biomass reaction found!")
    # Run pFBA
    solution = pfba(model)
    # Get growth rate from biomass reaction flux
    growth_rate = solution.fluxes[biomass_id]
    print("pFBA growth rate (biomass flux):", growth_rate)
    # Save all fluxes
    with open("pfba_fluxes.txt", "w") as f:
        for r in model.reactions:
            f.write(f"{r.id} : {solution.fluxes[r.id]}\n")
    print("pFBA fluxes saved to pfba_fluxes.txt.")

if __name__ == "__main__":
    main()