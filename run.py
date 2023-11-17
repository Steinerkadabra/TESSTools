from tess_tools.lc_tractor.extract_shortcadence import download_SC
from tess_tools.lc_interact.plot_lc import plot_lc
from tess_tools.lc_interact.phasephold_lc import phasefold_savgol

# download_SC("TIC 98545929", [33, 34], folder= "star_data")

phasefold_savgol("star_data/RS Cha", init_period=1.675, period_change=0.001, savgol_window = 300)
plot_lc("star_data/RS Cha")