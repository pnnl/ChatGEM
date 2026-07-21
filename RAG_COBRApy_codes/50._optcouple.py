# Description:
# This script performs an OptCouple-based strain design analysis on the iJN1463 genome-scale
# metabolic model using COBRApy and the StrainDesign (straindesign) library. The purpose of
# OptCouple is to identify reaction knockouts that enforce a minimum level of growth–product
# coupling between biomass formation and succinate production, ensuring that succinate
# synthesis is obligatorily linked to cellular growth. Cameo package must not be used. DO NOT
# USE CAMEO, use SDmodule from this codebase only
#
# The script begins by loading the SBML-formatted metabolic model (iJN1463.xml) and applying
# environmental constraints. These include fixing the glucose uptake rate and restricting
# selected exchange reactions to secretion-only operation. In addition, a small set of
# internal reactions is completely blocked to simulate secretion-only or pathway shutdown
# conditions as defined by the user.
#
# An exclusion list of reactions is defined to prevent OptCouple from selecting essential
# maintenance and respiratory reactions as knockout candidates, ensuring biologically
# meaningful strain designs.
#
# The biomass reaction (BIOMASS_KT2440_WT3) is used as the inner objective representing cellular
# growth, while the succinate exchange reaction (EX_succ_e) is defined as the target product
# reaction. The script first computes the wild-type optimal growth rate and flux distribution.
# A separate optimization is then performed to calculate the theoretical maximum succinate
# production rate without imposing growth requirements.
#
# An OptCouple module is constructed with a specified minimum growth-coupling potential
# (min_gcp), which enforces that succinate production must scale with biomass formation by at
# least this coupling strength. This formulation identifies reaction deletions that guarantee
# product formation whenever the organism grows.
#
# The OptCouple problem is solved allowing up to one reaction knockout and returning only the
# best-ranked solution. Candidate knockout sets are screened to remove any solution containing
# reactions from the exclusion list. The first valid solution is selected and used to generate
# a mutant metabolic model.
#
# The mutant strain is constructed by knocking out the selected reactions and is simulated
# under the same environmental conditions as the wild-type model to obtain mutant growth and
# flux distributions.
#
# The results are exported to an Excel file (OptCouple_iJN1463_succinate.xlsx) containing:
#  - A summary sheet reporting wild-type and mutant growth rates, succinate production rates,
#    theoretical maximum succinate yield, the number of knockouts, and the enforced minimum
#    growth-coupling potential (min_gcp).
#  - A knockouts sheet listing all reactions selected for deletion by OptCouple.
#  - A flux_comparison sheet containing reaction-wise wild-type and mutant fluxes and their
#    differences (mutant minus wild type).
#
# A helper function is included to correctly interpret knockout indicators returned by the
# StrainDesign solver. The script automatically selects the first valid OptCouple solution
# that satisfies exclusion constraints and coupling requirements.
#
# Sample ChatGPT prompts to get this code:
# - Python script for OptCouple strain design using COBRApy and StrainDesign.
# - How to enforce growth-coupled succinate production in iJN1463 using OptCouple.
# - Example of OptCouple-based reaction knockout analysis in Python.
# - COBRApy + StrainDesign workflow for growth–product coupling optimization.
# - Export OptCouple mutant vs wild-type fluxes to Excel using pandas.

import cobra
import straindesign as sd
import pandas as pd

def is_knockout_value(val):
    if isinstance(val, str):
        return val.strip().upper() in {"KO", "KNOCKOUT", "DEL", "DELETE"}
    try:
        return val < 0
    except Exception:
        return False

def main():
    # ----------------------
    # Load model
    # ----------------------
    model = cobra.io.read_sbml_model("iJN1463.xml")

    # ----------------------
    # Set glucose uptake
    # ----------------------
    if "EX_glc__D_e" in model.reactions:
        model.reactions.get_by_id("EX_glc__D_e").lower_bound = -6.0

    # ----------------------
    # Set other exchanges to secretion-only
    # ----------------------
    secretion_only = [
        "EX_co2_e",
        "EX_glcn_e",
        "EX_25dkglcn_e",
        "EX_catechol_e",
        "EX_lac__D_e",
        "EX_meoh_e",
        "EX_4oxptn_e",
        "EX_ac_e"
    ]
    for rxn_id in secretion_only:
        if rxn_id in model.reactions:
            r = model.reactions.get_by_id(rxn_id)
            r.lower_bound = 0.0
            r.upper_bound = 1000.0

    # ----------------------
    # Block secretion ONLY (commented out as per your code)
    # ----------------------
    block_only = ["FUM",
                  "ATPS4rpp",
                  "HCO3E",
                  "H2CO3D"
    ]
    for rxn_id in block_only:
        if rxn_id in model.reactions:
            r = model.reactions.get_by_id(rxn_id)
            r.lower_bound = 0.0
            r.upper_bound = 0.0

    # ----------------------
    # EXCLUSION LIST - prevent OptKnock from touching these
    # ----------------------
    excluded_rxns = [
        "ATPM",
        "CAT",
        "CYO1_KT",
        "CYTBO3_4pp",
        "CYTCAA3pp",
        "NADH16pp",
        "NAt3_1p5pp",
        "PItex",
        "PPK"
    ]

    # ----------------------
    # Biomass & target
    # ----------------------
    biomass_rxn = "BIOMASS_KT2440_WT3"
    #biomass_rxn = "BIOMASS_Ecoli_core_w_GAM"
    target_rxn = "EX_succ_e"

    # verify target exists
    if target_rxn not in model.reactions:
        raise ValueError(f"Exchange reaction {target_rxn} not found in the model.")

    # ----------------------
    # Wild-type baseline
    # ----------------------
    model.objective = biomass_rxn
    wt_sol = model.optimize()
    wt_flux = wt_sol.fluxes
    wt_growth = wt_sol.objective_value

    # ----------------------
    # Theoretical max succinate (unconstrained)
    # ----------------------
    succ_model = model.copy()
    succ_model.objective = target_rxn
    theoretical_succ = succ_model.optimize().objective_value

    # ----------------------
    # OptCouple module definition
    # ----------------------
    min_gcp_value = 0.01  # coupling strength requirement (tune as needed)
    optcouple_module = sd.SDModule(
        model,
        sd.OPTCOUPLE,
        inner_objective=biomass_rxn,
        prod_id=target_rxn,
        min_gcp=min_gcp_value
    )

    # ----------------------
    # Run OptCouple
    # ----------------------
    max_solutions = 1
    solutions = sd.compute_strain_designs(
        model,
        sd_modules=[optcouple_module],
        max_solutions=max_solutions,    # number of solutions solutions
        max_cost=1,                     # number of knockouts knockouts
        solution_approach=sd.BEST
    )

    # ----------------------
    # Loop over solutions, pick first one not containing excluded KOs
    # ----------------------
    ko_list = []
    found_valid = False
    for i in range(max_solutions):
        try:
            reaction_sd = solutions.get_reaction_sd(i)
        except Exception:
            break  # No more solutions
        candidate_ko = []
        skip = False
        for rd in reaction_sd:
            for rxn_id, val in rd.items():
                if is_knockout_value(val):
                    if rxn_id in excluded_rxns:
                        print(f"Solution {i}: Skipped because knockout {rxn_id} is in exclusion list.")
                        skip = True
                    candidate_ko.append(rxn_id)
        if skip:
            continue  # Try next solution
        ko_list = candidate_ko
        found_valid = True
        print(f"Selected solution {i} with knockouts: {ko_list}")
        break  # Only first valid solution

    if not found_valid or not ko_list:
        print("No valid OptCouple solution found (all contain excluded reactions).")
        return

    # ----------------------
    # Build mutant from solution
    # ----------------------
    mutant = model.copy()
    for rxn_id in ko_list:
        if rxn_id in mutant.reactions:
            mutant.reactions.get_by_id(rxn_id).knock_out()

    # ----------------------
    # Mutant simulation
    # ----------------------
    mutant.objective = biomass_rxn
    mut_sol = mutant.optimize()
    mut_flux = mut_sol.fluxes

    # ----------------------
    # Prepare Excel output
    # ----------------------
    wt_df = wt_flux.reset_index()
    wt_df.columns = ["reaction_id", "wt_flux"]
    mut_df = mut_flux.reset_index()
    mut_df.columns = ["reaction_id", "mutant_flux"]
    compare_df = wt_df.merge(mut_df, on="reaction_id", how="outer")
    compare_df["delta(mut-wt)"] = compare_df["mutant_flux"] - compare_df["wt_flux"]

    ko_df = pd.DataFrame({"knockout_reaction_id": ko_list})

    summary_df = pd.DataFrame({
        "metric": [
            "WT_growth",
            "WT_succinate",
            "Mut_growth",
            "Mut_succinate",
            "theoretical_max_succinate",
            "num_KOs",
            "min_gcp"
        ],
        "value": [
            wt_growth,
            wt_flux.get(target_rxn, float("nan")),
            mut_sol.objective_value,
            mut_flux.get(target_rxn, float("nan")),
            theoretical_succ,
            len(ko_list),
            min_gcp_value
        ]
    })

    out_file = "OptCouple_iJN1463_succinate.xlsx"
    with pd.ExcelWriter(out_file, engine="openpyxl") as writer:
        summary_df.to_excel(writer, sheet_name="summary", index=False)
        ko_df.to_excel(writer, sheet_name="knockouts", index=False)
        compare_df.to_excel(writer, sheet_name="flux_comparison", index=False)

    print(f"OptCouple computation complete. Results written to: {out_file}")

if __name__ == "__main__":
    main()