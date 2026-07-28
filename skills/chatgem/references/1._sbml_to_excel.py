# Description:
# This script loads a genome-scale metabolic model from SBML (XML) format using COBRApy
# and converts its structure into a multi-sheet Excel file for downstream analysis.
# The output Excel contains sheets for model info, reactions, metabolites, genes,
# compartments, and (optionally) the stoichiometric matrix.
# This workflow helps bridge model data with more accessible spreadsheet analysis.
#
# Sample ChatGPT prompts to get this code:
# - Python code to convert SBML XML metabolic model to multi-sheet Excel with COBRApy and pandas.
# - Export SBML model contents to Excel (reactions/metabolites/genes/compartments/matrix) in Python.
# - COBRApy script to read XML and create annotated Excel workbook.
# - How do I save SBML model stoichiometry and details to Excel in Python?
# - Convert a metabolic model from XML to Excel (multi-tab) using COBRApy and pandas.

import os
import cobra
import pandas as pd
from cobra.io import read_sbml_model
_ASSETS = os.path.join(os.path.dirname(__file__), '..', 'assets')

def convert_xml_to_excel(xml_file_path, excel_file_path):
    """
    Convert SBML XML model to Excel file with multiple sheets
    Parameters:
    xml_file_path (str): Path to the XML file
    excel_file_path (str): Path for output Excel file
    """
    # Load the model
    print("Loading model from XML file...")
    model = read_sbml_model(xml_file_path)
    print(f"Model loaded successfully!")
    print(f"Model ID: {model.id}")
    print(f"Number of reactions: {len(model.reactions)}")
    print(f"Number of metabolites: {len(model.metabolites)}")
    print(f"Number of genes: {len(model.genes)}")
    # Create Excel writer object
    with pd.ExcelWriter(excel_file_path, engine='openpyxl') as writer:
        # 1. Model Information Sheet
        model_info = {
            'Property': ['Model ID', 'Model Name', 'Number of Reactions',
                        'Number of Metabolites', 'Number of Genes', 'Objective'],
            'Value': [model.id, model.name, len(model.reactions),
                     len(model.metabolites), len(model.genes),
                     str(model.objective)]
        }
        pd.DataFrame(model_info).to_excel(writer, sheet_name='Model_Info', index=False)
        # 2. Reactions Sheet
        print("Processing reactions...")
        reactions_data = []
        for reaction in model.reactions:
            reactions_data.append({
                'ID': reaction.id,
                'Name': reaction.name,
                'Equation': reaction.build_reaction_string(),
                'Lower_Bound': reaction.lower_bound,
                'Upper_Bound': reaction.upper_bound,
                'Objective_Coefficient': reaction.objective_coefficient,
                'Subsystem': reaction.subsystem,
                'Gene_Rule': reaction.gene_reaction_rule,
                'Notes': str(reaction.notes) if reaction.notes else '',
                'Annotation': str(reaction.annotation) if reaction.annotation else ''
            })
        reactions_df = pd.DataFrame(reactions_data)
        reactions_df.to_excel(writer, sheet_name='Reactions', index=False)
        # 3. Metabolites Sheet
        print("Processing metabolites...")
        metabolites_data = []
        for metabolite in model.metabolites:
            metabolites_data.append({
                'ID': metabolite.id,
                'Name': metabolite.name,
                'Formula': metabolite.formula,
                'Charge': metabolite.charge,
                'Compartment': metabolite.compartment,
                'Notes': str(metabolite.notes) if metabolite.notes else '',
                'Annotation': str(metabolite.annotation) if metabolite.annotation else ''
            })
        metabolites_df = pd.DataFrame(metabolites_data)
        metabolites_df.to_excel(writer, sheet_name='Metabolites', index=False)
        # 4. Genes Sheet
        print("Processing genes...")
        genes_data = []
        for gene in model.genes:
            # Get reactions associated with this gene
            associated_reactions = [r.id for r in gene.reactions]
            genes_data.append({
                'ID': gene.id,
                'Name': gene.name,
                'Associated_Reactions': ', '.join(associated_reactions),
                'Number_of_Reactions': len(associated_reactions),
                'Notes': str(gene.notes) if gene.notes else '',
                'Annotation': str(gene.annotation) if gene.annotation else ''
            })
        genes_df = pd.DataFrame(genes_data)
        genes_df.to_excel(writer, sheet_name='Genes', index=False)
        # 5. Compartments Sheet
        print("Processing compartments...")
        compartments_data = []
        for comp_id, comp_name in model.compartments.items():
            # Count metabolites in each compartment
            met_count = len([m for m in model.metabolites if m.compartment == comp_id])
            compartments_data.append({
                'ID': comp_id,
                'Name': comp_name,
                'Number_of_Metabolites': met_count
            })
        compartments_df = pd.DataFrame(compartments_data)
        compartments_df.to_excel(writer, sheet_name='Compartments', index=False)
        # 6. Stoichiometric Matrix (optional - might be large)
        print("Creating stoichiometric matrix...")
        try:
            # Get stoichiometric matrix
            S = cobra.util.create_stoichiometric_matrix(model)
            S_df = pd.DataFrame(
                S,  # Don't use .todense()!
                index=[m.id for m in model.metabolites],
                columns=[r.id for r in model.reactions]
            )
            S_df.to_excel(writer, sheet_name='Stoichiometric_Matrix')
        except Exception as e:
            print(f"Warning: Could not create stoichiometric matrix sheet: {e}")
    print(f"Excel file saved successfully as: {excel_file_path}")

# Main execution
if __name__ == "__main__":
    # File paths
    xml_file = os.path.join(_ASSETS, 'iAF1260.xml')  # Input XML file
    excel_file = "iAF1260.xlsx"  # Output Excel file
    try:
        convert_xml_to_excel(xml_file, excel_file)
        print("Conversion completed successfully!")
    except FileNotFoundError:
        print(f"Error: Could not find the file '{xml_file}'. Please check the file path.")
    except Exception as e:
        print(f"Error during conversion: {e}")