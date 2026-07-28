# Description:
# This script loads the iAF1260 E. coli metabolic model using COBRApy, then removes a specific
# metabolite (identified by BiGG ID, here "glc__D_c") from the model object.
# This operation deletes the metabolite and all associated reactions containing it,
# which is useful for model curation, editing, or simulating metabolite deletions.

# Sample ChatGPT prompts to get this code:
# - Python code to remove a metabolite by ID from a COBRApy metabolic model.
# - COBRApy: how can I delete a metabolite from my model object?
# - Example for deleting a specific metabolite in SBML model using COBRApy.
# - Remove a metabolite and associated reactions by BiGG ID in COBRApy and print confirmation.
# - COBRApy script to remove and check a metabolite in a loaded metabolic model.

import os
import cobra
_ASSETS = os.path.join(os.path.dirname(__file__), '..', 'assets')
# Load model
model = cobra.io.read_sbml_model(os.path.join(_ASSETS, 'iAF1260.xml'))
# Remove metabolite by BiGG ID
model.metabolites.get_by_id("glc__D_c").remove_from_model()
print("Metabolite removed")