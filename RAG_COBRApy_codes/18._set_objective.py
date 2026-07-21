# Description:
# This script loads a metabolic model from an SBML file ("iAF1260.xml") using COBRApy,
# sets the objective function to a specific biomass reaction (replace `biomass_id` as needed),
# and prints the new objective expression. This is essential for defining the metabolic goal
# in simulations, like setting biomass production as the model objective.

# Sample ChatGPT prompts to get this code:
# - Python code to set the objective to a biomass reaction in a COBRApy model.
# - COBRApy: how can I change the objective to a different reaction?
# - Set and print the objective function for a metabolic model using COBRApy.
# - COBRApy script to update the objective to a specific reaction.
# - Change model objective and show updated expression in Python COBRApy.

import cobra
def main():
    model = cobra.io.read_sbml_model("iAF1260.xml")
    # Replace with the actual biomass reaction ID you found
    biomass_id = "BIOMASS_Ec_iAF1260_core_59p81M"
    model.objective = biomass_id
    print("Objective set to:", model.objective.expression)

if __name__ == "__main__":
    main()