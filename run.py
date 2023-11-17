from tess_tools.lc_tractor.extract_shortcadence import download_SC
from tess_tools.lc_interact.plot_lc import plot_lc

download_SC("TIC 98545929", [33, 34], folder= "star_data")

# plot_lc("star_data/RS Cha")