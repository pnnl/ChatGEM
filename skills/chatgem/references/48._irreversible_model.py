import os
import cobra
from cobra.io import read_sbml_model, write_sbml_model
from cobra import Reaction
import pandas as pd
import numpy as np
_ASSETS = os.path.join(os.path.dirname(__file__), '..', 'assets')

def is_exchange_reaction(rxn):
    """Check if a reaction is an exchange reaction"""
    # Check by ID patterns
    if any(prefix in rxn.id.upper() for prefix in ['EX_', 'DM_', 'SINK_', 'DEMAND_']):
        return True
    # Check if reaction has only one metabolite (typical for exchange reactions)
    if len(rxn.metabolites) == 1:
        return True
    # Check if all metabolites are in extracellular compartment
    if all('_e' in met.id or '_e0' in met.id or 'extracellular' in met.compartment.lower()
          for met in rxn.metabolites):
        return True
    return False

def split_reversible_reactions(model):
    model_copy = model.copy()
    reactions_to_remove = []
    reactions_to_add = []
    
    # Different bounds for regular vs exchange reactions
    regular_ub = 1000.0
    exchange_ub = 10.0
    default_lb = 0.0
    
    for rxn in model_copy.reactions:
        # Determine appropriate upper bound
        ub = exchange_ub if is_exchange_reaction(rxn) else regular_ub
        
        if rxn.lower_bound < 0 and rxn.upper_bound > 0:
            # Forward reaction
            fwd = Reaction(f"{rxn.id}_FWD")
            fwd.name = f"{rxn.name} (fwd)"
            fwd.lower_bound, fwd.upper_bound = default_lb, ub
            fwd.gene_reaction_rule = rxn.gene_reaction_rule
            fwd.subsystem = rxn.subsystem
            fwd.add_metabolites(rxn.metabolites)
            
            # Reverse reaction
            rev = Reaction(f"{rxn.id}_REV")
            rev.name = f"{rxn.name} (rev)"
            rev.lower_bound, rev.upper_bound = default_lb, ub
            rev.gene_reaction_rule = rxn.gene_reaction_rule
            rev.subsystem = rxn.subsystem
            rev.add_metabolites({met: -coeff for met, coeff in rxn.metabolites.items()})
            
            reactions_to_add.extend([fwd, rev])
            reactions_to_remove.append(rxn)
            
        elif rxn.lower_bound < 0 and rxn.upper_bound <= 0:
            # Convert reverse-only to forward
            new_rxn = Reaction(f"{rxn.id}_FWD")
            new_rxn.name = f"{rxn.name} (conv)"
            new_rxn.lower_bound, new_rxn.upper_bound = default_lb, ub
            new_rxn.gene_reaction_rule = rxn.gene_reaction_rule
            new_rxn.subsystem = rxn.subsystem
            new_rxn.add_metabolites({met: -coeff for met, coeff in rxn.metabolites.items()})
            
            reactions_to_add.append(new_rxn)
            reactions_to_remove.append(rxn)
            
        elif rxn.lower_bound >= 0 and rxn.upper_bound > 0:
            # Already irreversible forward reaction - just update bounds
            rxn.lower_bound = default_lb
            rxn.upper_bound = ub
    
    model_copy.remove_reactions(reactions_to_remove)
    model_copy.add_reactions(reactions_to_add)
    
    # Final pass to ensure all reactions have correct bounds
    for rxn in model_copy.reactions:
        if rxn.lower_bound < 0:
            rxn.lower_bound = default_lb
        # Set appropriate upper bound based on reaction type
        target_ub = exchange_ub if is_exchange_reaction(rxn) else regular_ub
        if rxn.upper_bound != target_ub:
            rxn.upper_bound = target_ub
    
    return model_copy

def run_custom_optimization():
    # Load model
    model = read_sbml_model(os.path.join(_ASSETS, 'iJB1325.xml'))
    
    # Split reversible reactions
    split_model = split_reversible_reactions(model)
    
    # Save irreversible model
    write_sbml_model(split_model, "irreviJB1325.xml")
    print("Irreversible model saved as: irreviJB1325.xml")
    
    # Find GROWTH reaction
    growth_reaction_id = None
    growth_reaction = None
    if "GROWTH" in [r.id for r in split_model.reactions]:
        growth_reaction_id = "GROWTH"
        growth_reaction = split_model.reactions.get_by_id("GROWTH")
    elif "GROWTH_FWD" in [r.id for r in split_model.reactions]:
        growth_reaction_id = "GROWTH_FWD"
        growth_reaction = split_model.reactions.get_by_id("GROWTH_FWD")
    else:
        # Look for any reaction containing "GROWTH" in the name
        for rxn in split_model.reactions:
            if "GROWTH" in rxn.id.upper():
                growth_reaction_id = rxn.id
                growth_reaction = rxn
                break
    
    if not growth_reaction:
        print("GROWTH reaction not found!")
        return None, None, None
    
    # Set up custom objective: max GROWTH - 0.0001 * (sum of all reactions)
    # First clear any existing objective
    split_model.objective = {}
    
    # Create objective dictionary using reaction objects
    objective_dict = {}
    
    for rxn in split_model.reactions:
        if rxn.id == growth_reaction_id:
            # GROWTH: maximize (1.0) - minimize as part of sum (0.0001)
            objective_dict[rxn] = 1.0 - 0.0001
        else:
            # All other reactions: minimize as part of sum
            objective_dict[rxn] = -0.0001
    
    # Set the objective using reaction objects
    split_model.objective = objective_dict
    
    # Optimize
    solution = split_model.optimize()
    
    if solution.status == 'optimal':
        growth_flux = solution.fluxes[growth_reaction_id]
        objective_value = solution.objective_value
        
        # Calculate total flux sum using numpy
        total_flux_sum = np.sum(np.abs(solution.fluxes.values))
        
        print(f"Optimization successful!")
        print(f"GROWTH reaction ({growth_reaction_id}) flux: {growth_flux:.6f}")
        print(f"Objective value: {objective_value:.6f}")
        print(f"Total flux sum: {total_flux_sum:.6f}")
        
        # Create Excel data for flux distribution
        data = []
        for r in split_model.reactions:
            is_exchange = is_exchange_reaction(r)
            flux_value = solution.fluxes[r.id]
            # Get objective coefficient
            obj_coeff = 1.0 - 0.0001 if r.id == growth_reaction_id else -0.0001
            
            data.append({
                'Reaction_ID': r.id,
                'Name': r.name,
                'Subsystem': r.subsystem,
                'Genes': r.gene_reaction_rule,
                'Equation': r.reaction,
                'Lower_Bound': r.lower_bound,
                'Upper_Bound': r.upper_bound,
                'Is_Exchange': is_exchange,
                'Flux': float(flux_value),  # Convert to Python float
                'Abs_Flux': float(np.abs(flux_value)),  # Use numpy abs and convert to float
                'Active': abs(flux_value) > 1e-6,
                'Objective_Coefficient': obj_coeff
            })
        
        df = pd.DataFrame(data)
        
        # Save Excel file
        with pd.ExcelWriter("flux_iCH384.xlsx") as writer:
            # All fluxes
            df.to_excel(writer, sheet_name='All_Fluxes', index=False)
            
            # Active reactions only
            active_df = df[df['Active']]
            active_df.to_excel(writer, sheet_name='Active_Only', index=False)
            
            # Exchange reactions
            exchange_df = df[df['Is_Exchange'] == True]
            if not exchange_df.empty:
                exchange_df.to_excel(writer, sheet_name='Exchange_Reactions', index=False)
            
            # Regular reactions
            regular_df = df[df['Is_Exchange'] == False]
            if not regular_df.empty:
                regular_df.to_excel(writer, sheet_name='Regular_Reactions', index=False)
            
            # Summary
            summary_data = {
                'Metric': [
                    'GROWTH_Flux', 
                    'Objective_Value', 
                    'Total_Flux_Sum',
                    'Total_Reactions', 
                    'Active_Reactions',
                    'Exchange_Reactions', 
                    'Regular_Reactions'
                ],
                'Value': [
                    float(growth_flux),
                    float(objective_value),
                    float(total_flux_sum),
                    len(split_model.reactions),
                    int(df['Active'].sum()),
                    len(exchange_df),
                    len(regular_df)
                ]
            }
            pd.DataFrame(summary_data).to_excel(writer, sheet_name='Summary', index=False)
        
        print("Flux distribution saved to: flux_iJB1325.xlsx")
        
        return growth_flux, objective_value, total_flux_sum
    else:
        print(f"Optimization failed with status: {solution.status}")
        return None, None, None

# Run the analysis
growth_flux, objective_value, total_flux_sum = run_custom_optimization()