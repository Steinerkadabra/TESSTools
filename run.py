from tess_tools.lc_tractor.extract_shortcadence import download_SC
from tess_tools.lc_tractor.extract_longcadence import get_lc_from_FFI
from tess_tools.lc_interact.plot_lc import plot_lc
from tess_tools.lc_interact.phasephold_lc import phasefold_savgol


from astroquery.gaia import Gaia

download_SC("TIC 179040251", [61], folder= "star_data")

# phasefold_savgol("star_data/RS Cha", init_period=1.675, period_change=0.001, savgol_window = 300)
# plot_lc("star_data/RS Cha")


# get_lc_from_FFI("TIC 40607177",  "40607177", [61], folder= "star_data/ngc_1901", maglim=3)
# get_lc_from_FFI("TIC 40607177",  "40607177", [61, 62, 63, 64, 65, 66, 67, 68, 69], folder= "star_data/ngc_1901", maglim=3)
