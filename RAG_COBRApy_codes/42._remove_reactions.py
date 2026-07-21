# Description:
# This script loads the iAF1260 E. coli metabolic model using COBRApy, then removes a specific 
# reaction (by BiGG reaction ID, here "PGI") from the model. This operation permanently deletes
# the reaction from the loaded model object, which can be useful for simulating reaction knockouts
# or manually curating and editing genome-scale metabolic models.

# Sample ChatGPT prompts to get this code:
# - Python code to remove a reaction by ID from a COBRApy metabolic model.
# - COBRApy: how do I delete a reaction from the model object?
# - Example of removing a specific reaction from SBML model using COBRApy.
# - Delete a reaction by BiGG ID in COBRApy and print confirmation.
# - COBRApy script to remove and check a reaction in the loaded model.

import cobra
# Load model
model = cobra.io.read_sbml_model("iAF1260.xml")
# Remove reaction by BiGG ID
model.reactions.get_by_id("PGI").remove_from_model()
print("Reaction removed")