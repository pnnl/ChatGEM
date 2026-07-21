import cobra
import pandas as pd
import numpy as np
import gurobipy as gp
from gurobipy import GRB

# ============================================================
# PARAMETERS
# ============================================================
PARAMS = {
    "MODEL_SBML": "iJN1463_irrev.xml",
    "PROTEIN_XLSX": "MW_kcat.xlsx",
    "UPPER_BOUND_XLSX": "upper_bound.xlsx",

    "BIOMASS_RXN": "BIOMASS_KT2440_WT3",
    "OUTPUT_XLSX": "ecGEM_output.xlsx",

    "PROTEOME_UB_g_per_gDW": 0.20,
    "PROTEIN_COST_DIVISOR": 1000.0,
    "FLUX_PENALTY": 0.00,

    "gurobi_OutputFlag": 0,
    "gurobi_DualReductions": 0,
}

# ============================================================
# MODEL BUILD
# ============================================================
def build_and_solve(model, prot_cost, upper_bounds, params):
    m = gp.Model("ecGEM")
    m.setParam("OutputFlag", params["gurobi_OutputFlag"])
    m.setParam("DualReductions", params["gurobi_DualReductions"])

    # --------------------------------------------------------
    # VARIABLES
    # --------------------------------------------------------
    v = {}
    for rxn in model.reactions:
        if rxn.id in upper_bounds:
            ub = max(0.0, upper_bounds[rxn.id])
        else:
            ub = max(0.0, float(rxn.upper_bound))
        v[rxn.id] = m.addVar(lb=0, ub=ub, name=f"v_{rxn.id}")
    m.update()

    # --------------------------------------------------------
    # MASS BALANCE
    # --------------------------------------------------------
    for met in model.metabolites:
        m.addConstr(
            gp.quicksum(float(rxn.metabolites[met]) * v[rxn.id] for rxn in met.reactions) == 0,
            name=f"mass_{met.id}",
        )

    # --------------------------------------------------------
    # PROTEOME CONSTRAINT  (<=, not ==)
    # --------------------------------------------------------
    if prot_cost:
        prot_expr = gp.quicksum(float(cost) * v[rid] for rid, cost in prot_cost.items()) / params["PROTEIN_COST_DIVISOR"]
        m.addConstr(prot_expr <= params["PROTEOME_UB_g_per_gDW"], name="protein_capacity")

    # --------------------------------------------------------
    # OBJECTIVE: maximize biomass - 0.001 * sum(all fluxes)
    # --------------------------------------------------------
    biomass_id = params["BIOMASS_RXN"]
    flux_penalty = params["FLUX_PENALTY"] * gp.quicksum(v[rid] for rid in v)
    m.setObjective(v[biomass_id] - flux_penalty, GRB.MAXIMIZE)
    m.optimize()

    return m, v

# ============================================================
# LOAD MODEL
# ============================================================
model = cobra.io.read_sbml_model(PARAMS["MODEL_SBML"])

if PARAMS["BIOMASS_RXN"] not in [r.id for r in model.reactions]:
    raise ValueError("Biomass reaction not found.")

# ============================================================
# LOAD UPPER BOUNDS FROM EXCEL
# ============================================================
ub_raw = pd.read_excel(PARAMS["UPPER_BOUND_XLSX"], usecols=[0, 1], names=["Reactions", "Upper_bound"])
model_rxns = set(r.id for r in model.reactions)
upper_bounds = {}

for _, row in ub_raw.iterrows():
    rid = str(row["Reactions"])
    if rid in model_rxns and pd.notna(row["Upper_bound"]):
        upper_bounds[rid] = float(row["Upper_bound"])

print(f"Loaded {len(upper_bounds)} reaction upper bounds from Excel.")

# ============================================================
# LOAD PROTEIN COST
# ============================================================
prot_raw = pd.read_excel(PARAMS["PROTEIN_XLSX"], usecols=[0, 1], names=["Reactions", "Protein_cost"])
prot_cost = {}

for _, row in prot_raw.iterrows():
    rid = str(row["Reactions"])          # FIX: was row["Protein_cost"] — wrong column
    if rid in model_rxns and pd.notna(row["Protein_cost"]):
        prot_cost[rid] = float(row["Protein_cost"])

print(f"Loaded {len(prot_cost)} protein costs from Excel.")

# ============================================================
# SOLVE
# ============================================================
m, v = build_and_solve(model, prot_cost, upper_bounds, PARAMS)

if m.status != GRB.OPTIMAL:
    raise RuntimeError("Optimization did not converge.")

# ============================================================
# EXTRACT SOLUTION
# ============================================================
fluxes = {rid: var.X for rid, var in v.items()}

# ============================================================
# DATAFRAMES
# ============================================================
flux_df = pd.DataFrame.from_dict(fluxes, orient="index", columns=["flux"]).reset_index()
flux_df.columns = ["Reaction", "flux"]

# ============================================================
# PROTEIN USAGE
# ============================================================
protein_rows = []
p_total = 0.0
for rid, cost in prot_cost.items():
    pr = cost * fluxes.get(rid, 0) / PARAMS["PROTEIN_COST_DIVISOR"]
    protein_rows.append({
        "Reaction": rid,
        "Protein_Required (g/gDW)": pr
    })
    p_total += pr

protein_df = pd.DataFrame(protein_rows)

# ============================================================
# WRITE EXCEL
# ============================================================
with pd.ExcelWriter(PARAMS["OUTPUT_XLSX"], engine="openpyxl") as writer:
    flux_df.to_excel(writer, sheet_name="flux", index=False)
    protein_df.to_excel(writer, sheet_name="protein_requirements", index=False)

print("Biomass growth rate (/hr) =", v[PARAMS["BIOMASS_RXN"]].X)
print(f"Total protein requirement (g/gDW): {p_total}")