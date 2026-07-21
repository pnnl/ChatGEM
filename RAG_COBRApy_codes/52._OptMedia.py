import cobra
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from optlang import Constraint, Variable, Objective

# ---------------------------
# PARAMETERS
# ---------------------------
MODEL_FILE = "iJN1463.xml"
BIOMASS_RXN_ID = "BIOMASS_KT2440_WT3"
PRODUCT_RXN_ID = "EX_succ_e"
BIOMASS_TARGET_FRAC = 0.5
PRODUCT_TARGET_FRAC = 0.5
BIG_M = 1000
MEDIA_UPTAKE_LB = -6.0
MEDIA_UPPER_BOUND = 1000.0
MINIMAL_MEDIA_FILE = "minimal_media.xlsx"
FLUX_DISTRIBUTION_FILE = "flux_distribution.xlsx"
GROWTH_FRACTION_FOR_CHECK = 1.0
TOL = 1e-9
VERBOSE = True

# ---------------------------
# FORCE MEDIA COMPONENTS
# ---------------------------
# Must always be present (available for uptake)
REQUIRED_MEDIA = ["EX_glc__D_e", "EX_nh4_e", "EX_so4_e"]
# Must never be used (uptake disallowed)
FORBIDDEN_MEDIA = ["EX_icit_e"]

# ---------------------------
# MEDIA EXCHANGES (grouped)
# ---------------------------
MEDIA_EXCHANGES = [
    # Carbon sources / central metabolites
    "EX_cell4_e","EX_cit_e","EX_fru_e","EX_fum_e","EX_glc__D_e","EX_glyc_e",
    "EX_glycol_e","EX_lac__L_e","EX_man_e","EX_meoh_e","EX_mal__D_e","EX_pyr_e",
    "EX_rib__D_e","EX_etoh_e","EX_indole_e","EX_skm_e","EX_akg_e","EX_icit_e",
    "EX_oaa_e",
    # Nitrogen sources and amino acids
    "EX_nh4_e","EX_no3_e","EX_urea_e","EX_orn_e","EX_ala_B_e","EX_ala__L_e",
    "EX_arg__L_e","EX_asn__L_e","EX_asp__L_e","EX_cys__L_e","EX_gln__L_e",
    "EX_glu__L_e","EX_gly_e","EX_his__L_e","EX_ile__L_e","EX_leu__L_e",
    "EX_lys__L_e","EX_met__L_e","EX_phe__L_e","EX_pro__L_e","EX_ser__L_e",
    "EX_thr__L_e","EX_trp__L_e","EX_tyr__L_e","EX_val__L_e",
    # Inorganic nutrients
    "EX_so4_e","EX_pi_e",
    # Metals / ions
    "EX_ca2_e","EX_cd2_e","EX_cobalt2_e","EX_cu2_e","EX_cu_e","EX_fe2_e",
    "EX_fe3_e","EX_mg2_e","EX_mn2_e","EX_mobd_e","EX_na1_e","EX_ni2_e",
    "EX_zn2_e","EX_k_e",
    # Other essentials / atmosphere
    "EX_h2o_e","EX_h_e","EX_o2_e"
]

# ---------------------------
# PLOTTING DEFAULTS
# ---------------------------
plt.rcParams.update({
    "font.size": 14,
    "axes.labelweight": "bold",
    "axes.linewidth": 2,
    "xtick.major.width": 2,
    "ytick.major.width": 2,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "savefig.bbox": "tight",
    "figure.dpi": 150,
})

def pnas_axes(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(width=2, length=6, direction="out")

# ---------------------------
# HELPERS
# ---------------------------
def log(msg: str):
    if VERBOSE:
        print(msg)

def apply_medium_from_y(m, y_sol):
    """
    Apply a binary medium selection y to a COBRA model:
      - All candidate exchanges: no uptake (LB=0).
      - Selected exchanges (y>0.5): uptake allowed down to MEDIA_UPTAKE_LB.
      - REQUIRED_MEDIA exchanges are ALWAYS treated as selected.
      - FORBIDDEN_MEDIA exchanges are NEVER allowed (forced closed).
    """
    required_set = set(REQUIRED_MEDIA)
    forbidden_set = set(FORBIDDEN_MEDIA)

    # Close uptake for all candidate exchanges
    for rid in MEDIA_EXCHANGES:
        if rid in m.reactions:
            rxn = m.reactions.get_by_id(rid)
            rxn.lower_bound = 0.0
            rxn.upper_bound = MEDIA_UPPER_BOUND

    # Open uptake for anything selected OR required, except forbidden
    for rid in MEDIA_EXCHANGES:
        if rid not in m.reactions:
            continue
        if rid in forbidden_set:
            continue
        is_required = rid in required_set
        is_selected = (y_sol.get(rid, 0.0) > 0.5)
        if is_required or is_selected:
            rxn = m.reactions.get_by_id(rid)
            rxn.lower_bound = MEDIA_UPTAKE_LB
            rxn.upper_bound = MEDIA_UPPER_BOUND

def robust_check(y_sol, alpha, growth_frac=1.0, tol=1e-9):
    """
    Robustness test:
      1) Maximize biomass => mu*
      2) Fix biomass at growth_frac*mu* and minimize product exchange => v_succ_min
      3) Robust if v_succ_min >= alpha*(growth_frac*mu*) - tol
    """
    m = cobra.io.read_sbml_model(MODEL_FILE)
    bio = m.reactions.get_by_id(BIOMASS_RXN_ID)
    prd = m.reactions.get_by_id(PRODUCT_RXN_ID)

    apply_medium_from_y(m, y_sol)

    # Step 1: maximize growth
    m.objective = bio
    sol1 = m.optimize()
    if sol1.status != "optimal":
        return False, None, None
    mu = sol1.objective_value

    # Step 2: fix growth and minimize product
    with m:
        bio.lower_bound = growth_frac * mu
        m.objective = prd
        m.objective_direction = "min"
        sol2 = m.optimize()
        if sol2.status != "optimal":
            return False, mu, None
        min_succ = sol2.objective_value

    # Step 3: check robustness condition
    ok = (min_succ >= alpha * (growth_frac * mu) - tol)
    return ok, mu, min_succ

# ---------------------------
# SANITY CHECKS
# ---------------------------
missing_req = set(REQUIRED_MEDIA) - set(MEDIA_EXCHANGES)
missing_forb = set(FORBIDDEN_MEDIA) - set(MEDIA_EXCHANGES)
overlap = set(REQUIRED_MEDIA) & set(FORBIDDEN_MEDIA)
if missing_req:
    raise ValueError(f"REQUIRED_MEDIA not in MEDIA_EXCHANGES: {sorted(missing_req)}")
if missing_forb:
    raise ValueError(f"FORBIDDEN_MEDIA not in MEDIA_EXCHANGES: {sorted(missing_forb)}")
if overlap:
    raise ValueError(f"Media cannot be both required and forbidden: {sorted(overlap)}")

# ---------------------------
# 1) Compute alpha from WT benchmarks
# ---------------------------
base = cobra.io.read_sbml_model(MODEL_FILE)
base.objective = base.reactions.get_by_id(BIOMASS_RXN_ID)
wt_growth = base.optimize().objective_value

base.objective = base.reactions.get_by_id(PRODUCT_RXN_ID)
max_prod = base.optimize().objective_value

alpha = PRODUCT_TARGET_FRAC * max_prod / wt_growth

log(f"WT biomass (max): {wt_growth:.6g}")
log(f"Max succinate (EX): {max_prod:.6g}")
log(f"Biomass target (lb): {BIOMASS_TARGET_FRAC * wt_growth:.6g}")
log(f"Coupling alpha: {alpha:.6g}")

# ---------------------------
# 2) Build MILP
# ---------------------------
model = cobra.io.read_sbml_model(MODEL_FILE)
biomass_rxn = model.reactions.get_by_id(BIOMASS_RXN_ID)
prod_rxn = model.reactions.get_by_id(PRODUCT_RXN_ID)

# Start by allowing all candidate exchanges to be uptaken (binary constraints will shut them off)
for rxn_id in MEDIA_EXCHANGES:
    if rxn_id in model.reactions:
        rxn = model.reactions.get_by_id(rxn_id)
        rxn.lower_bound = MEDIA_UPTAKE_LB
        rxn.upper_bound = MEDIA_UPPER_BOUND

# Enforce biomass lower bound
model.add_cons_vars([Constraint(biomass_rxn.flux_expression, lb=BIOMASS_TARGET_FRAC * wt_growth)])

# Enforce coupling: v_product - alpha * v_biomass >= 0
prod_coupling = prod_rxn.flux_expression - alpha * biomass_rxn.flux_expression
model.add_cons_vars([Constraint(prod_coupling, lb=0)])

binary_vars = {}
required_set = set(REQUIRED_MEDIA)
forbidden_set = set(FORBIDDEN_MEDIA)

for rxn_id in MEDIA_EXCHANGES:
    if rxn_id in model.reactions:
        rxn = model.reactions.get_by_id(rxn_id)

        # Binary var: y=1 means nutrient is available for uptake
        y = Variable(f"y_{rxn_id}", type="binary")
        binary_vars[rxn_id] = y

        # If y=0 => v >= 0 (no uptake). If y=1 => v can be negative down to -BIG_M
        model.add_cons_vars([Constraint(rxn.flux_expression + BIG_M * y, lb=0)])

        # Force required nutrients to always be present (y == 1)
        if rxn_id in required_set:
            model.add_cons_vars([Constraint(y, lb=1, ub=1)])

        # Force forbidden nutrients to never be present (y == 0)
        if rxn_id in forbidden_set:
            model.add_cons_vars([Constraint(y, lb=0, ub=0)])

# Objective: minimize number of nutrients in the medium (sum of binaries)
model.objective = Objective(sum(binary_vars.values()), direction="min")

# ---------------------------
# 3) Cutting-plane loop
# ---------------------------
max_iters = 500
y_sol = None
ok = False

iterations_fig = []
mu_fig = []
min_succ_fig = []
margin_fig = []
media_count_fig = []
media_matrix_fig = pd.DataFrame(0, index=MEDIA_EXCHANGES, columns=[])

for it in range(1, max_iters + 1):
    status = model.solver.optimize()
    if status != "optimal":
        raise RuntimeError("MILP infeasible or solver failed.")

    # Read solution (binary selection)
    y_sol = {rid: model.solver.variables[var.name].primal for rid, var in binary_vars.items()}

    # Robustness check
    ok, mu, min_succ = robust_check(y_sol, alpha, growth_frac=GROWTH_FRACTION_FOR_CHECK, tol=TOL)

    iterations_fig.append(it)
    mu_fig.append(mu if mu is not None else 0.0)
    min_succ_fig.append(min_succ if min_succ is not None else 0.0)

    # Robustness margin
    if (mu is not None) and (min_succ is not None):
        margin = min_succ - alpha * (GROWTH_FRACTION_FOR_CHECK * mu)
    else:
        margin = 0.0
    margin_fig.append(margin)

    # Media selection vector for plotting (required are always 1; forbidden always 0)
    media_vec = []
    for rid in MEDIA_EXCHANGES:
        if rid in forbidden_set:
            media_vec.append(0)
        elif rid in required_set:
            media_vec.append(1)
        else:
            media_vec.append(int(y_sol.get(rid, 0) > 0.5))
    media_matrix_fig[it] = media_vec
    media_count_fig.append(int(sum(media_vec)))

    log(
        f"Iter {it:03d} | mu*={mu:.6g} | v_succ_min={min_succ:.6g} "
        f"| margin={margin:.6g} | robust={ok} | media_n={media_count_fig[-1]}"
    )

    if ok:
        break

    # Add a no-good cut to exclude this exact binary pattern next time
    on = [binary_vars[rid] for rid, v in y_sol.items() if v > 0.5]
    off = [binary_vars[rid] for rid, v in y_sol.items() if v <= 0.5]
    model.add_cons_vars(Constraint(sum(on) - sum(off), ub=len(on) - 1))

# ---------------------------
# 4) Save selected medium
# ---------------------------
# Ensure required nutrients are included; ensure forbidden are excluded
selected = sorted((set(REQUIRED_MEDIA) | {rid for rid, v in y_sol.items() if v > 0.5}) - set(FORBIDDEN_MEDIA))
pd.DataFrame({"Reaction": selected}).to_excel(MINIMAL_MEDIA_FILE, index=False)
log(f"Saved selected medium to: {MINIMAL_MEDIA_FILE} (n={len(selected)})")

# ---------------------------
# 5) Final FBA and flux export
# ---------------------------
final_model = cobra.io.read_sbml_model(MODEL_FILE)
apply_medium_from_y(final_model, y_sol)

final_model.objective = final_model.reactions.get_by_id(BIOMASS_RXN_ID)
fba_solution = final_model.optimize()
fba_solution.fluxes.to_excel(FLUX_DISTRIBUTION_FILE)
log(f"Saved flux distribution to: {FLUX_DISTRIBUTION_FILE}")

# ---------------------------
# 6) Figures (no titles; TIFF only)
# ---------------------------
# Figure 1
fig, ax1 = plt.subplots(figsize=(9, 5))
ax1.plot(iterations_fig, mu_fig, linewidth=2.8, color="black", label="Growth rate")
ax1.set_xlabel("Iterations")
ax1.set_ylabel("Growth rate (h$^{-1}$)")
pnas_axes(ax1)

ax2 = ax1.twinx()
ax2.plot(iterations_fig, min_succ_fig, linestyle="--", linewidth=2.8, color="#2b8cbe", label="Min succinate flux")
ax2.set_ylabel("EX_succ_e (mmol gDW$^{-1}$ h$^{-1}$)")
ax2.tick_params(axis="y", width=2, length=6, direction="out")
ax2.spines["top"].set_visible(False)

lines, labels = [], []
for ax in (ax1, ax2):
    h, l = ax.get_legend_handles_labels()
    lines += h
    labels += l

ax1.legend(
    lines, labels,
    loc="lower center",
    bbox_to_anchor=(0.5, 1.02),
    ncol=2,
    frameon=False,
    fontsize=12,
    handlelength=3
)

plt.tight_layout()
plt.savefig("fig1.tiff", dpi=300)
plt.show()

# Figure 2
fig, ax = plt.subplots(figsize=(8.5, 5))
ax.plot(iterations_fig, margin_fig, color="black", linewidth=2.0, zorder=1)
sc = ax.scatter(
    iterations_fig, margin_fig,
    c=media_count_fig, cmap="viridis",
    s=90, edgecolor="black", linewidth=0.8, zorder=2
)
ax.axhline(0, color="black", linewidth=2)
ax.set_xlabel("Iterations")
ax.set_ylabel("Robustness margin")
pnas_axes(ax)

cbar = plt.colorbar(sc, ax=ax, pad=0.02)
cbar.set_label("Media components (#)")
cbar.outline.set_linewidth(2)

plt.tight_layout()
plt.savefig("fig2.tiff", dpi=300)
plt.show()

# Figure 3: Media selection heatmap over iterations (VERTICAL) — no empty y labels
binary_matrix = media_matrix_fig.T  # rows=iterations, cols=exchanges
binary_matrix = binary_matrix.loc[:, (binary_matrix.sum(axis=0) > 0)]  # keep only ever-selected exchanges
binary_matrix_v = binary_matrix.T  # rows=components, cols=iterations
binary_matrix_v.index = [c.replace("EX_", "").replace("_e", "") for c in binary_matrix_v.index]

fig, ax = plt.subplots(
    figsize=(max(6, binary_matrix_v.shape[1] * 0.35),
             max(6, binary_matrix_v.shape[0] * 0.25))
)
sns.heatmap(
    binary_matrix_v,
    cmap=["#d73027", "#1a9850"],
    cbar=False,
    linewidths=0.6,
    linecolor="black",
    ax=ax
)
ax.set_xlabel("Iterations")
ax.tick_params(width=2, length=8, direction="out")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# Force all y ticks/labels (prevents blanks due to auto-subsampling)
ax.set_yticks(range(binary_matrix_v.shape[0]))
ax.set_yticklabels(binary_matrix_v.index, rotation=0)

# Optional: force all x ticks too
ax.set_xticks(range(binary_matrix_v.shape[1]))
ax.set_xticklabels(binary_matrix_v.columns, rotation=0)

plt.tight_layout()
plt.savefig("fig3.tiff", dpi=300)
plt.show()