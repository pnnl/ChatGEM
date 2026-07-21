import cobra
import pandas as pd
from optlang import Objective

# Load the model
model = cobra.io.read_sbml_model("irreviJB1325.xml")
print("Model loaded successfully")

# Load protein cost data
an_data = pd.read_excel("per_unit_protein_data_30C.xlsx", usecols=[0, 1], names=["Reactions", "Protein_cost"])
model_reaction_ids = set(r.id for r in model.reactions)
an_protein_costs = {}
for _, row in an_data.iterrows():
    reaction_id = row["Reactions"]
    if reaction_id in model_reaction_ids:
        an_protein_costs[reaction_id] = row["Protein_cost"]
print("Protein costs loaded")

# Validate GROWTH reaction
if "GROWTH" in model.reactions:
    growth_reaction = model.reactions.get_by_id("GROWTH")
else:
    print("GROWTH reaction not found!")
    exit()

print("Setting metabolic proteome constraint (metabolic protein only)...")

# Set up metabolic protein constraints
an_variables = {}
for reaction_id, cost in an_protein_costs.items():
    reaction = model.reactions.get_by_id(reaction_id)
    an_variables[reaction.forward_variable] = cost

# Total metabolic protein expression (from reaction costs)
# Divide by 1000 to get units g/gDW — constraint is still 600!
met_expr = sum(coeff * var for var, coeff in an_variables.items()) / 1000
proteome_constraint = model.problem.Constraint(
    met_expr, ub=0.60, name="AN_protein_capacity"
)
model.add_cons_vars(proteome_constraint)
print("Metabolic protein capacity constraint applied")

# Set objective: maximize GROWTH - 0.0001*(sum of all fluxes)
print("Setting objective: GROWTH - 0.0001*(sum of all fluxes)")
all_flux_expr = sum(reaction.forward_variable for reaction in model.reactions)
objective_expr = growth_reaction.forward_variable - 0.0001 * all_flux_expr
model.objective = Objective(objective_expr, direction='max')

# Optimize
print("Optimizing...")
solution = model.optimize()
if solution.status == "optimal":
    print(f"Optimization successful!")
    print(f"Growth rate: {growth_reaction.flux:.5f} per hour")
    print(f"Objective value: {solution.objective_value:.5f}")

    # === Everything reported in g/gDW ===
    an_growth = growth_reaction.flux
    p_met = sum(model.reactions.get_by_id(rxn).flux * cost for rxn, cost in an_protein_costs.items()) / 1000   # [g/gDW]
    total_flux_sum = sum(abs(reaction.flux) for reaction in model.reactions)

    print("\n-- PROTEIN REQUIREMENT (AN/an) --")
    print(f"Growth rate: {an_growth:.3f} per hour")
    print(f"Metabolic protein requirement (= total): {p_met:.3f} g/gDW")
    print(f"Total used: {p_met:.3f} g/gDW (should be ~0.6 g/gDW)")
    print(f"Sum of all fluxes: {total_flux_sum:.3f}")
else:
    print(f"Optimization failed with status: {solution.status}")
    exit()

# Export to Excel — everything in g/gDW ONLY
print("Exporting results...")
with pd.ExcelWriter("AN_Results_30C.xlsx") as writer:
    # Flux distribution
    an_flux_data = []
    for reaction_id in an_protein_costs.keys():
        if reaction_id in model.reactions:
            flux = model.reactions.get_by_id(reaction_id).flux
            an_flux_data.append({
                "Reaction_ID": reaction_id,
                "Flux (mmol/gDW/hr)": flux
            })
    an_flux_df = pd.DataFrame(an_flux_data)
    an_flux_df.to_excel(writer, sheet_name='Flux_Distribution', index=False)

    # Protein requirements by reaction (all /1000, so g/gDW)
    an_protein_data = []
    for reaction_id, cost in an_protein_costs.items():
        flux = model.reactions.get_by_id(reaction_id).flux
        protein_required = (cost * flux) / 1000  # g/gDW
        an_protein_data.append({
            "Reaction_ID": reaction_id,
            "Protein_Required (g/gDW)": protein_required
        })
    # Add metabolic total (= total)
    an_protein_data.append({
        "Reaction_ID": "TOTAL (metabolic only)",
        "Protein_Required (g/gDW)": p_met
    })
    an_protein_df = pd.DataFrame(an_protein_data)
    an_protein_df.to_excel(writer, sheet_name='Protein_Requirements', index=False)

    # Protein pool summary (ONLY g/gDW)
    summary_df = pd.DataFrame([
        {"Pool": "Metabolic (=Total)", "g/gDW": p_met}
    ])
    summary_df.to_excel(writer, sheet_name='Protein_Pools', index=False)

    # Overall summary
    overall_summary_df = pd.DataFrame([
        {"Metric": "Growth Rate (per hour)", "Value": an_growth},
        {"Metric": "Objective Value", "Value": solution.objective_value},
        {"Metric": "Total Protein Usage (g/gDW)", "Value": p_met},
        {"Metric": "Protein Capacity (g/gDW, for interpretation)", "Value": 0.6},
        {"Metric": "Total Flux Sum", "Value": total_flux_sum}
    ])
    overall_summary_df.to_excel(writer, sheet_name='Summary', index=False)

print("Results exported to AN_Results_30C.xlsx")