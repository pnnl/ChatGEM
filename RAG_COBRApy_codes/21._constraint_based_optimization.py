# Description:
# This script loads a metabolic model from an SBML file ("iAF1260.xml") using COBRApy,
# sets new bounds for the glucose exchange reaction (EX_glc__D_e), specifically allowing
# uptake up to 20 units and no secretion, then optimizes the model with these constraints.
# The resulting new biomass/growth rate is printed. This workflow demonstrates how to modify
# nutrient constraints and observe effects on growth in genome-scale models.

# Sample ChatGPT prompts to get this code:
# - Python code to constrain glucose uptake in a COBRApy model and print growth rate.
# - How to change exchange reaction bounds and optimize with COBRApy (example)?
# - COBRApy: set reaction bounds and get new objective value after optimization.
# - Example of modifying uptake bounds then running FBA using COBRApy.
# - Modify substrate constraints and show impact on COBRApy model growth.

import cobra
def main():
    # Load the model
    model = cobra.io.read_sbml_model("iAF1260.xml")
    # Set new bounds on glucose exchange reaction
    reaction = model.reactions.get_by_id("EX_glc__D_e")
    reaction.lower_bound = -20  # Allow uptake up to 20 units
    reaction.upper_bound = 0    # Prevent secretion
    # Optimize the model
    solution = model.optimize()
    # Print the new growth rate (objective value)
    print("New growth rate after constraint:", solution.objective_value)

if __name__ == "__main__":
    main()