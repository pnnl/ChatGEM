# Description:
# This script loads a COBRApy metabolic model (e_coli_core.xml), optionally allows model manipulation,
# and then exports the model as SBML to a temporary file. It runs Memote to perform a comprehensive
# quality and consistency check, saving the results as a HTML report. This workflow is useful for
# benchmarking, curating, and improving genome-scale metabolic models.
# Sample ChatGPT prompts to get this code:
# - Python code to run Memote quality report on a COBRApy SBML model and save results to HTML.
# - How do I benchmark a metabolic model with Memote using Python and COBRApy?
# - Example for exporting model as SBML and running Memote suite report.
# - COBRApy & Memote: automated report generation for model quality.
# - Generate Memote test report for genome-scale metabolic models in Python.

import cobra
import memote.suite.api
from pathlib import Path

model = cobra.io.read_sbml_model("e_coli_core.xml")

# Set the solver to GLPK exact
model.solver = 'glpk_exact'

# You can manipulate the model here if desired

# Save to a temporary SBML fashion (optional, but not necessary for Memote API)
cobra.io.write_sbml_model(model, "temp_model.xml")

# Run Memote: pass the model object!
results = memote.suite.api.test_model(model)

# Write the HTML report
memote.suite.api.snapshot_report(results, filename="memote_report.html")

print("Memote report generated: memote_report.html")