# Description:
# This script performs dynamic flux balance analysis (dFBA) using COBRApy and SciPy's ODE integrator on the iAF1260 E. coli model.
# It simulates growth on glucose over time, updating the glucose uptake rate with Michaelis-Menten kinetics in each ODE step.
# Results (biomass and glucose concentrations over time) are plotted and saved to a CSV file.
# This workflow demonstrates how to couple genome-scale metabolic models with dynamic environments for growth simulations.

# Sample ChatGPT prompts to get this code:
# - Python code for dynamic FBA simulation with COBRApy and glucose uptake kinetics.
# - COBRApy dynamic simulation example (dFBA, time-dependent growth and substrate depletion).
# - How do I combine ODE integration with FBA for E. coli growth in COBRApy?
# - Example script for Michaelis-Menten constrained dFBA using COBRApy and SciPy.
# - Dynamic FBA: simulate substrate-limited batch growth with COBRApy and save results.

import os
import numpy as np
from tqdm import tqdm
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import cobra
import pandas as pd

_ASSETS = os.path.join(os.path.dirname(__file__), "..", "assets")

# === 1. Load your genome-scale model ===
model = cobra.io.read_sbml_model(os.path.join(_ASSETS, "iAF1260.xml"))

# === 2. Identify relevant reaction IDs ===
# Biomass reaction (find one automatically)
biomass_rxn_id = None
for r in model.reactions:
    if "BIOMASS" in r.id.upper():
        biomass_rxn_id = r.id
        print("Biomass reaction:", biomass_rxn_id)
        break
if not biomass_rxn_id:
    raise ValueError("No biomass reaction found!")
# Glucose exchange reaction
glc_rxn_id = "EX_glc__D_e"


# === 3. Dynamic bounds function ===
def add_dynamic_bounds(model, y):
    """Use external concentrations to bound the uptake flux of glucose."""
    biomass, glucose = y
    # Michaelis-Menten kinetic uptake (change parameters as needed)
    vmax = 10  # maximum uptake rate (mmol/gDW/hr)
    km = 5  # Michaelis constant (mmol/L)
    glucose_max_import = -vmax * glucose / (km + glucose)
    model.reactions.get_by_id(glc_rxn_id).lower_bound = glucose_max_import


# === 4. Dynamic system function ===
def dynamic_system(t, y):
    """Calculate the time derivatives of biomass and glucose."""
    biomass, glucose = y
    with model:
        add_dynamic_bounds(model, y)
        solution = model.optimize()
        mu = solution.fluxes[biomass_rxn_id]  # 1/hr
        glc_uptake = -solution.fluxes[glc_rxn_id]  # mmol/gDW/hr
    # Rate of change of biomass and glucose
    d_biomass = mu * biomass
    d_glucose = -glc_uptake * biomass
    # Progress bar (optional for notebook)
    if dynamic_system.pbar is not None:
        dynamic_system.pbar.update(1)
        dynamic_system.pbar.set_description("t = {:.3f}".format(t))
    return [d_biomass, d_glucose]


dynamic_system.pbar = None


# === 5. Feasibility event function ===
def infeasible_event(t, y):
    """Detect simulation infeasibility or glucose exhaustion."""
    biomass, glucose = y
    if glucose <= 0.01:  # stop if glucose nearly gone
        return 0
    with model:
        add_dynamic_bounds(model, y)
        solution = model.optimize()
        mu = solution.fluxes[biomass_rxn_id]
    return mu - infeasible_event.epsilon


infeasible_event.epsilon = 1e-6
infeasible_event.direction = 1
infeasible_event.terminal = True

# === 6. Run the dynamic FBA simulation ===
ts = np.linspace(0, 10, 101)  # Simulate for 10 hours
y0 = [0.1, 10]  # Initial biomass (gDW/L), glucose (mmol/L)
with tqdm(total=len(ts)) as pbar:
    dynamic_system.pbar = pbar
    sol = solve_ivp(
        fun=dynamic_system,
        events=[infeasible_event],
        t_span=(ts.min(), ts.max()),
        y0=y0,
        t_eval=ts,
        rtol=1e-6,
        atol=1e-8,
        method="BDF",
    )

# === 7. Plotting results ===
plt.figure(figsize=(8, 4))
plt.plot(sol.t, sol.y[0], label="Biomass (gDW/L)")
plt.plot(sol.t, sol.y[1], label="Glucose (mmol/L)")
plt.xlabel("Time (h)")
plt.ylabel("Concentration")
plt.title("Dynamic FBA - iAF1260 E. coli Growth on Glucose")
plt.legend()
plt.tight_layout()
plt.show()

# === 8. Save results to file ===
df = pd.DataFrame(
    {"Time (h)": sol.t, "Biomass (gDW/L)": sol.y[0], "Glucose (mmol/L)": sol.y[1]}
)
df.to_csv("dfba_iaf1260_results.txt", sep="\t", index=False)
print("Results saved to dfba_iaf1260_results.txt")
