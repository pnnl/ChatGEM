import pandas as pd
import re
import os
import matplotlib.pyplot as plt
from cobra import io
from cobra.exceptions import OptimizationError
import matplotlib.patches as mpatches

# File paths for SBML model, Excel output, and transcriptomics data
xml_filepath = os.path.join(os.getcwd(), 'e_coli_core.xml')
excel_filepath = os.path.join(os.getcwd(), 'Model.xlsx')
transcriptomics_filepath = os.path.join(os.getcwd(), 'transcriptomics_data.xlsx')

# Load the SBML model using COBRApy
model = io.read_sbml_model(xml_filepath)
try:
    # Optimize the model based on its default objective function to find the baseline solution
    solution = model.optimize()
except OptimizationError:
    # Catch optimization errors, if any, to ensure the program continues execution
    pass

# Collect reaction data including bounds and GPR rules
reaction_data = []
for reaction in model.reactions:
    lower_bound = reaction.lower_bound
    upper_bound = reaction.upper_bound
    gpr = str(reaction.gene_reaction_rule)  # Convert GPR rule to string
    reaction_data.append({
        'Reaction ID': reaction.id,
        'Lower Bound': lower_bound,
        'Upper Bound': upper_bound,
        'GPR': gpr
    })
reaction_df = pd.DataFrame(reaction_data)  # Convert reaction data to a DataFrame

# Extract unique genes from the model for gene expression matching
unique_genes = {gene.id for gene in model.genes}

# Load transcriptomics data from an Excel file
transcriptomics_data = pd.read_excel(transcriptomics_filepath)

# Convert the transcriptomics data into a dictionary for quick lookup of gene expression levels
expression_map = transcriptomics_data.set_index('Gene_ID').T.to_dict('list')
condition_names = transcriptomics_data.columns[1:]  # Extract condition names (columns excluding 'Gene_ID')

# Prepare gene expression data for model genes, filling in 0 for missing genes
gene_data = []
for gene_id in unique_genes:
    expression_values = expression_map.get(gene_id, [0] * len(condition_names))
    gene_data.append({'Gene_ID': gene_id, **dict(zip(condition_names, expression_values))})
gene_df = pd.DataFrame(gene_data)  # Convert gene expression data to a DataFrame
gene_df = gene_df[['Gene_ID'] + list(condition_names)]  # Ensure the column order includes condition names

def evaluate_gpr(gpr, expression_map, condition_index):
    """Evaluate the GPR rule for a reaction based on gene expression data for a given condition."""
    # Extract genes from the GPR string
    genes = re.findall(r"[\w]+", gpr)
    
    # If GPR is a single gene and expression data is available, return its expression value
    if len(genes) == 1 and genes[0] in expression_map:
        return expression_map[genes[0]][condition_index]
    
    # Default to 1 if there is no GPR rule
    if not gpr:
        return 1
    
    # Replace logical operators with Python equivalents
    gpr = gpr.replace('and', '&&').replace('or', '||')
    parts = re.split(r"\(|\)", gpr)  # Split expression by parentheses
    
    stack = []  # Stack to evaluate expressions based on logical operators
    for part in parts:
        sub_exprs = re.split(r"&&|\|\|", part)  # Split expression by logical operators
        
        values = []  # List to hold expression values
        for sub_expr in sub_exprs:
            sub_expr = sub_expr.strip()  # Clean whitespace
            if sub_expr in expression_map:
                values.append(expression_map[sub_expr][condition_index])  # Add expression value if available
            elif re.findall(r"[\w]+", sub_expr):
                values.append(0)  # Default expression value to 0 if gene not found
        
        # Evaluate expressions based on logical "and" and "or"
        if '&&' in part and values:
            stack.append(min(values))
        elif '||' in part and values:
            stack.append(max(values))
    
    return max(stack) if stack else 1  # Default to 1 if no part is valid

# Calculate preliminary bound scores based on GPR evaluations for each condition
preliminary_bound_scores = []
for reaction in model.reactions:
    if reaction.gene_reaction_rule:
        score_entry = {'Reaction ID': reaction.id, 'GPR': reaction.gene_reaction_rule}
        for cond_index, condition_name in enumerate(condition_names):
            score = evaluate_gpr(reaction.gene_reaction_rule, expression_map, cond_index)
            score_entry[condition_name] = score
        preliminary_bound_scores.append(score_entry)
preliminary_bound_df = pd.DataFrame(preliminary_bound_scores)  # Convert scores to DataFrame

# Normalize the preliminary bound data to scale GPR scores
normalized_preliminary_bound_df = preliminary_bound_df.copy()
for condition_name in condition_names:
    max_value = normalized_preliminary_bound_df[condition_name].max()
    if max_value != 0:
        normalized_preliminary_bound_df[condition_name] = (
            normalized_preliminary_bound_df[condition_name] / max_value * 1000
        )

def create_condition_tabs(condition_name, reactions_df, normalized_bound_df):
    """Create condition-specific reaction bounds, treating zero GPR score as if without GPR constraints."""
    condition_df = reactions_df[['Reaction ID', 'Lower Bound', 'Upper Bound']].copy()
    for idx, row in condition_df.iterrows():
        if row['Reaction ID'] in normalized_bound_df['Reaction ID'].values:
            normalized_value = normalized_bound_df.loc[normalized_bound_df['Reaction ID'] == row['Reaction ID'], condition_name].values[0]
            
            if normalized_value == 0:
                # If GPR score is 0, treat the reaction without a GPR, using original bounds
                condition_df.at[idx, 'Lower Bound'] = row['Lower Bound']
                condition_df.at[idx, 'Upper Bound'] = row['Upper Bound']
            else:
                # Scale bounds according to normalized GPR score
                if row['Lower Bound'] < 0:  # If the lower bound is negative
                    condition_df.at[idx, 'Lower Bound'] = -normalized_value
                condition_df.at[idx, 'Upper Bound'] = normalized_value

        # Handle special case for exchange reactions, assumed to have an 'EX_' prefix
        if "EX_" in row['Reaction ID']:
            if row['Lower Bound'] < 0:  # Reversible exchange reaction: set bounds
                condition_df.at[idx, 'Lower Bound'] = -10
                condition_df.at[idx, 'Upper Bound'] = 10
            else:  # Forward-only exchange reaction: set bounds
                condition_df.at[idx, 'Lower Bound'] = 0
                condition_df.at[idx, 'Upper Bound'] = 10
    
    return condition_df

# Output all data into an Excel file with multiple sheets to store various stages of calculated data
with pd.ExcelWriter(excel_filepath, engine='openpyxl') as writer:
    # Write essential tables and data to different sheets
    reaction_df.to_excel(writer, index=False, sheet_name='model_reactions')
    gene_df.to_excel(writer, index=False, sheet_name='expression')
    preliminary_bound_df.to_excel(writer, index=False, sheet_name='preliminary_bound')
    normalized_preliminary_bound_df.to_excel(writer, index=False, sheet_name='normalized_preliminary_bound')

    # Generate condition-specific sheets
    for condition_name in condition_names:
        condition_specific_df = create_condition_tabs(condition_name, reaction_df, normalized_preliminary_bound_df)
        condition_specific_df.to_excel(writer, index=False, sheet_name=condition_name)

# Read the Excel file for model evaluation under different conditions
xls = pd.ExcelFile(excel_filepath)

# Identify which sheets correspond to condition-specific reaction bounds
start_index = xls.sheet_names.index("normalized_preliminary_bound") + 1

# Track biomass growth rates for each condition
biomass_growth_rates = []

# Evaluate model under the condition-specific reaction bounds
for sheet_name in xls.sheet_names[start_index:]:
    df = pd.read_excel(xls, sheet_name=sheet_name)
    # Reload original SBML model for each condition evaluation step
    model = io.read_sbml_model(xml_filepath)

    # Update the model with condition-specific bounds
    for index, row in df.iterrows():
        reaction_id = row['Reaction ID']  # Obtain reaction ID from the Excel sheet
        lower_bound = row['Lower Bound']  # Read lower bound for the reaction
        upper_bound = row['Upper Bound']  # Read upper bound for the reaction

        # Implement updated bounds in the model
        if reaction_id in model.reactions:
            reaction = model.reactions.get_by_id(reaction_id)
            reaction.lower_bound = lower_bound
            reaction.upper_bound = upper_bound

    # Optimize and record biomass growth rate
    try:
        solution = model.optimize()
        biomass_growth_rates.append((sheet_name, solution.objective_value))
    except OptimizationError:
        # Record None if optimization fails
        biomass_growth_rates.append((sheet_name, None))

    # Save the adjusted model as a new SBML file named after the condition
    output_filepath = os.path.join(os.getcwd(), f"{sheet_name}.xml")
    io.write_sbml_model(model, output_filepath)
    print(f"Updated model saved as: {output_filepath}")

# Plot biomass growth rates across conditions
condition_names = [name for name, _ in biomass_growth_rates]
growth_rates = [rate if rate is not None else 0 for name, rate in biomass_growth_rates]  # Handle None values appropriately

# Setup color cycle for bar chart
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
colors *= int((len(biomass_growth_rates) / len(colors)) + 1)  # Repeat colors to ensure enough colors are available

plt.figure(figsize=(12, 8))
bars = plt.bar(condition_names, growth_rates, color=colors[:len(condition_names)], edgecolor='black')

# Configure plot appearance
plt.xlabel('Conditions', fontsize=12, fontweight='bold')
plt.ylabel('Biomass Growth Rate', fontsize=12, fontweight='bold')
plt.title('Biomass Growth Rate by Condition', fontsize=14, fontweight='bold')
plt.xticks(rotation=45, ha='right', fontsize=10, fontweight='bold')
plt.yticks(fontsize=10, fontweight='bold')
plt.tight_layout()

# Optional: Include legend to describe conditions
legend_patches = [mpatches.Patch(color=colors[i], label=condition_names[i]) for i in range(len(condition_names))]
plt.legend(handles=legend_patches, title='Conditions', fontsize=10, title_fontsize='12', loc='upper right')

# Display the plot
plt.show()