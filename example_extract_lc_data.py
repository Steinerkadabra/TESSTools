# This example shows how to extract a light curve for RS Cha in the sectors
# 64, 65, and 66 from the Full Frame Images  and save it in the folder 
# example_data/download_sc_data. The folder must exist beforehand. RS Cha
# has short cadence data available, but is taking for the porpuse of showing 
# the savgol filter in example_use_savgol_filter_to_remove_binary_signal.py

from tess_tools.lc_tractor.extract_longcadence import get_lc_from_FFI


get_lc_from_FFI("TIC 323292655",  "323292655", [64, 65, 66], folder= "example_data/extract_lc_data", maglim=3)
# The aperture needed for the correction to work properly here is huge. Keep in mind, the BCSAP flux here is fine.
# See example_use_savgol_filter_to_remove_binary_signal.py