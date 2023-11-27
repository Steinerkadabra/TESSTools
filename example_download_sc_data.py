# This example shows how to download short cadence data for Beta Pic in the sectors
# 32, 33, and 34 and save it in the folder example_data/download_sc_data. The folder 
# must exist beforehand.

from tess_tools.lc_tractor.extract_shortcadence import download_SC


download_SC("Beta Pic", [32, 33, 34], folder= "example_data/download_sc_data")