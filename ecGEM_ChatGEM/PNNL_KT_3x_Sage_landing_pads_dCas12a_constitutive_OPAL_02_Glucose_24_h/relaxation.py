import cobra
import pandas as pd
import numpy as np
import gurobipy as gp
from gurobipy import GRB

PARAMS = {
    "MODEL_SBML": "iJN1463_irrev.xml",
    "PROTEIN_XLSX": "MW_kcat.xlsx",
    "UPPER_BOUND_XLSX": "upper_bound.xlsx",
    "BIOMASS_RXN": "BIOMASS_KT2440_WT3",
    "OUTPUT_XLSX": "relaxation_reactions.xlsx",
    "PROTEOME_UB_g_per_gDW": 0.60,
    "PROTEIN_COST_DIVISOR": 1000.0,
    "FLUX_PENALTY": 0.001,
    "gurobi_OutputFlag": 0,
    "gurobi_DualReductions": 0,
}

def build_and_solve(model, prot_cost, upper_bounds, params):
    m = gp.Model("ecGEM")
    m.setParam("OutputFlag", params["gurobi_OutputFlag"])
    m.setParam("DualReductions", params["gurobi_DualReductions"])
    v = {}
    for rxn in model.reactions:
        if rxn.id in upper_bounds:
            ub = max(0.0, upper_bounds[rxn.id])
        else:
            ub = max(0.0, float(rxn.upper_bound))
        v[rxn.id] = m.addVar(lb=0, ub=ub, name=f"v_{rxn.id}")
    m.update()
    for met in model.metabolites:
        m.addConstr(
            gp.quicksum(float(rxn.metabolites[met]) * v[rxn.id] for rxn in met.reactions) == 0,
            name=f"mass_{met.id}",
        )
    if prot_cost:
        prot_expr = gp.quicksum(float(cost) * v[rid] for rid, cost in prot_cost.items()) / params["PROTEIN_COST_DIVISOR"]
        m.addConstr(prot_expr <= params["PROTEOME_UB_g_per_gDW"], name="protein_capacity")
    biomass_id = params["BIOMASS_RXN"]
    flux_penalty = params["FLUX_PENALTY"] * gp.quicksum(v[rid] for rid in v)
    m.setObjective(v[biomass_id] - flux_penalty, GRB.MAXIMIZE)
    m.optimize()
    return m, v

model = cobra.io.read_sbml_model(PARAMS["MODEL_SBML"])
if PARAMS["BIOMASS_RXN"] not in [r.id for r in model.reactions]:
    raise ValueError("Biomass reaction not found.")

model_rxns = set(r.id for r in model.reactions)

ub_raw = pd.read_excel(PARAMS["UPPER_BOUND_XLSX"], usecols=[0, 1], names=["Reactions", "Upper_bound"])
upper_bounds = {}
for _, row in ub_raw.iterrows():
    rid = str(row["Reactions"])
    if rid in model_rxns and pd.notna(row["Upper_bound"]):
        upper_bounds[rid] = float(row["Upper_bound"])
print(f"Loaded {len(upper_bounds)} reaction upper bounds.")

prot_raw = pd.read_excel(PARAMS["PROTEIN_XLSX"], usecols=[0, 1], names=["Reactions", "Protein_cost"])
prot_cost = {}
for _, row in prot_raw.iterrows():
    rid = str(row["Reactions"])
    if rid in model_rxns and pd.notna(row["Protein_cost"]):
        prot_cost[rid] = float(row["Protein_cost"])
print(f"Loaded {len(prot_cost)} protein costs.")

print("\nSolving unconstrained model...")
m_free, v_free = build_and_solve(model, prot_cost, upper_bounds={}, params=PARAMS)
if m_free.status != GRB.OPTIMAL:
    raise RuntimeError("Unconstrained model is not optimal.")
free_growth = v_free[PARAMS["BIOMASS_RXN"]].X
print(f"Unconstrained growth: {free_growth:.4f}")

bottlenecks = {}
for rid, ub in upper_bounds.items():
    free_flux = v_free[rid].X if rid in v_free else 0.0
    if ub > 0 and (free_flux - ub) / ub > 1.0:
        bottlenecks[rid] = {
            "Original_UB": ub,
            "Free_flux": free_flux,
            "Relative_deviation": (free_flux - ub) / ub
        }
print(f"Found {len(bottlenecks)} bottleneck reactions.")

bottleneck_df = pd.DataFrame([
    {
        "Reaction": rid,
        "Original_UB": vals["Original_UB"],
        "Free_flux": round(vals["Free_flux"], 6),
        "Relative_deviation": round(vals["Relative_deviation"], 4),
        "New_UB": 1000.0
    }
    for rid, vals in bottlenecks.items()
]).sort_values("Relative_deviation", ascending=False)

bottleneck_df.to_excel(PARAMS["OUTPUT_XLSX"], index=False)
print(f"Reactions relaxed: {len(bottlenecks)}")
print(f"Results written to {PARAMS['OUTPUT_XLSX']}")

print("\nSolving relaxed model...")
final_ubs = upper_bounds.copy()
for rid in bottlenecks:
    final_ubs[rid] = 1000.0
m, v = build_and_solve(model, prot_cost, final_ubs, PARAMS)
if m.status != GRB.OPTIMAL:
    raise RuntimeError("Relaxed model is not optimal.")
print(f"Unconstrained growth: {free_growth:.4f}")
print(f"Constrained growth (after relaxation): {v[PARAMS['BIOMASS_RXN']].X:.4f}")