# Description:
# This script loads a metabolic model from an SBML file ("iAF1260.xml") using COBRApy,
# runs Flux Balance Analysis (FBA), and saves the growth rate (objective value) and each reaction's flux
# to a text file called "fba_fluxes.txt". This helps document or analyze the simulation’s results.

# Sample ChatGPT prompts to get this code:
# - Python code to run FBA and save all fluxes from a COBRApy model.
# - COBRApy: how to print objective value and write fluxes to a file?
# - Save FBA results (growth rate + reaction fluxes) after optimization using COBRApy.
# - COBRApy example to export FBA solution to a text file.
# - Write FBA objective and fluxes for all reactions to file with COBRApy.

import cobra
def main():
    model = cobra.io.read_sbml_model("iAF1260.xml")
    solution = model.optimize()
    with open("fba_fluxes.txt", "w") as f:
        f.write(f"Objective value (growth rate): {solution.objective_value}\n")
        for r in model.reactions:
            f.write(f"{r.id} : {solution.fluxes[r.id]}\n")
    print("FBA results saved to fba_fluxes.txt.")

if __name__ == "__main__":
    main()