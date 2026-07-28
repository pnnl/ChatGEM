import os
import cobra
import pandas as pd
import numpy as np
import gurobipy as gp
from gurobipy import GRB
_ASSETS = os.path.join(os.path.dirname(__file__), '..', 'assets')

# ============================================================
# PARAMETERS
# ============================================================
PARAMS = {
    "MODEL_SBML": os.path.join(_ASSETS, 'iJN1463_irrev.xml'),
    "DELG_XLSX": os.path.join(_ASSETS, 'delG_post_opt.xlsx'),
    "PROTEIN_XLSX": os.path.join(_ASSETS, 'per_unit_protein_data.xlsx'),

    "BIOMASS_RXN": "BIOMASS_KT2440_WT3",
    "OUTPUT_XLSX": "ecOptMDFPathway.xlsx",

    "R_kJ_per_mol_K": 8.314e-3,
    "T_K": 298.15,

    "Cmin": 1e-6,
    "Cmax": 100.0,

    "PROTEOME_UB_g_per_gDW": 0.60,
    "PROTEIN_COST_DIVISOR": 1000.0,

    "lambda_B": 1e-4,
    "lambda_F": 1e-4,

    "min_biomass": 1e-6,
    "min_flux": 1e-6,

    "BIG_M": 1000,

    "gurobi_OutputFlag": 0,
    "gurobi_DualReductions": 0,
}

# ============================================================
# Gibbs calculation
# ============================================================
def compute_delG_from_solution(model, delG0_dict, lnC_vals, RT):
    delG = {}
    for rxn in model.reactions:
        s = 0.0
        for met, nu in rxn.metabolites.items():
            s += float(nu) * float(lnC_vals[met.id])
        delG[rxn.id] = float(delG0_dict[rxn.id]) + RT * s
    return delG

# ============================================================
# MODEL BUILD
# ============================================================
def build_and_solve(model, delG0, prot_cost, params):
    Rg = params["R_kJ_per_mol_K"]
    T = params["T_K"]
    RT = Rg * T

    x_lb = np.log(params["Cmin"])
    x_ub = np.log(params["Cmax"])
    BIG_M = params["BIG_M"]

    m = gp.Model("ecOptMDFPathway")
    m.setParam("OutputFlag", params["gurobi_OutputFlag"])
    m.setParam("DualReductions", params["gurobi_DualReductions"])

    # --------------------------------------------------------
    # VARIABLES
    # --------------------------------------------------------
    v = {}
    y = {}
    for rxn in model.reactions:
        ub = max(0.0, float(rxn.upper_bound))
        v[rxn.id] = m.addVar(lb=0, ub=ub, name=f"v_{rxn.id}")
        y[rxn.id] = m.addVar(vtype=GRB.BINARY, name=f"y_{rxn.id}")

    x = {met.id: m.addVar(lb=x_lb, ub=x_ub, name=f"x_{met.id}") for met in model.metabolites}
    f = {rxn.id: m.addVar(lb=-GRB.INFINITY, name=f"f_{rxn.id}") for rxn in model.reactions}
    B = m.addVar(lb=0, name="B")
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
    # THERMODYNAMICS
    # --------------------------------------------------------
    for rxn in model.reactions:
        dG0 = float(delG0[rxn.id])
        m.addConstr(
            f[rxn.id] == -(dG0 + RT * gp.quicksum(float(coeff) * x[met.id] for met, coeff in rxn.metabolites.items())),
            name=f"thermo_{rxn.id}",
        )

    # --------------------------------------------------------
    # FLUX-ACTIVITY COUPLING
    # --------------------------------------------------------
    for rxn in model.reactions:
        ub = max(0.0, float(rxn.upper_bound))
        m.addConstr(v[rxn.id] <= ub * y[rxn.id], name=f"flux_ub_{rxn.id}")
        m.addConstr(v[rxn.id] >= params["min_flux"] * y[rxn.id], name=f"flux_lb_{rxn.id}")

    # --------------------------------------------------------
    # MDF COUPLING
    # --------------------------------------------------------
    for rxn in model.reactions:
        m.addConstr(f[rxn.id] >= B - BIG_M * (1 - y[rxn.id]), name=f"mdf_{rxn.id}")

    # --------------------------------------------------------
    # PROTEOME CONSTRAINT
    # --------------------------------------------------------
    if prot_cost:
        prot_expr = gp.quicksum(float(cost) * v[rid] for rid, cost in prot_cost.items()) / params["PROTEIN_COST_DIVISOR"]
        m.addConstr(prot_expr <= params["PROTEOME_UB_g_per_gDW"], name="protein_capacity")

    # --------------------------------------------------------
    # MINIMUM BIOMASS
    # --------------------------------------------------------
    biomass_id = params["BIOMASS_RXN"]
    if params["min_biomass"] > 0:
        m.addConstr(v[biomass_id] >= params["min_biomass"], name="min_biomass")

    # --------------------------------------------------------
    # OBJECTIVE
    # --------------------------------------------------------
    total_flux = gp.quicksum(v[rxn.id] for rxn in model.reactions)
    m.setObjective(v[biomass_id] + params["lambda_B"] * B - params["lambda_F"] * total_flux, GRB.MAXIMIZE)
    m.optimize()

    return m, v, x, B, RT

# ============================================================
# LOAD MODEL
# ============================================================
model = cobra.io.read_sbml_model(PARAMS["MODEL_SBML"])

if PARAMS["BIOMASS_RXN"] not in [r.id for r in model.reactions]:
    raise ValueError("Biomass reaction not found.")

# ============================================================
# LOAD ΔG°
# ============================================================
delG_df = pd.read_excel(PARAMS["DELG_XLSX"])
delG0 = pd.Series(delG_df["delG0"].values, index=delG_df["Reaction"]).to_dict()

# ============================================================
# LOAD PROTEIN COST
# ============================================================
prot_raw = pd.read_excel(PARAMS["PROTEIN_XLSX"], usecols=[0, 1], names=["Reactions", "Protein_cost"])
model_rxns = set(r.id for r in model.reactions)
prot_cost = {}

for _, row in prot_raw.iterrows():
    rid = str(row["Reactions"])
    if rid in model_rxns and pd.notna(row["Protein_cost"]):
        prot_cost[rid] = float(row["Protein_cost"])

# ============================================================
# SOLVE
# ============================================================
m, v, x, B, RT = build_and_solve(model, delG0, prot_cost, PARAMS)

if m.status != GRB.OPTIMAL:
    raise RuntimeError("Optimization did not converge.")

# ============================================================
# EXTRACT SOLUTION
# ============================================================
fluxes = {rid: var.X for rid, var in v.items()}
lnC = {mid: var.X for mid, var in x.items()}
C = {mid: np.exp(val) for mid, val in lnC.items()}
delG_model = compute_delG_from_solution(model, delG0, lnC, RT)

# ============================================================
# DATAFRAMES
# ============================================================
flux_df = pd.DataFrame.from_dict(fluxes, orient="index", columns=["flux"]).reset_index()
flux_df.columns = ["Reaction", "flux"]

conc_df = pd.DataFrame.from_dict(C, orient="index", columns=["concentration"]).reset_index()
conc_df.columns = ["Metabolite", "concentration"]

thermo_df = pd.DataFrame({
    "Reaction": [r.id for r in model.reactions],
    "delG0_final": [delG0[r.id] for r in model.reactions],
    "delG_model": [delG_model[r.id] for r in model.reactions],
})

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
    conc_df.to_excel(writer, sheet_name="concentration", index=False)
    thermo_df.to_excel(writer, sheet_name="thermodynamics", index=False)
    protein_df.to_excel(writer, sheet_name="protein_requirements", index=False)

print("Biomass growth rate (/hr) =", v[PARAMS["BIOMASS_RXN"]].X)
print(f"Total protein requirement (g/gDW): {p_total}")
print("Minimum thermodynamic driving force (Kj/mol) =", B.X)