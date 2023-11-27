# This example showcases how to plot a previously downloaded or extracted
# light curve as well as how to use the phasefold_savgol function to remove
# the periodic signal of the ecplising binary signal in the BCSAP light curve
# of RS Cha. 

# Run example_extract_lc_data.py before running this example.

from tess_tools.lc_interact.phasephold_lc import phasefold_savgol
from tess_tools.lc_interact.plot_lc import plot_lc

# We call plot_lc with the folder to extraction results for RS Cha. We can first
# look as the comnbined corrected data. We will see, that the correction did not 
# work properly (If you have taken a huge aperture, it might have).

# # Write the number in the terminal that is before TIC 323292655_combined_corrected64_65_66.txt 
plot_lc("example_data/extract_lc_data/TIC 323292655")

# # Now, we call plot_lc again to look at the backgorund correct simple aperture photometry flux

# # Write the number in the terminal that is before TIC 323292655_combined_BCSAP64_65_66.txt 
plot_lc("example_data/extract_lc_data/TIC 323292655")

# # As you can see, the BCSAP is more then enough. In order to look at he pulsation light curve
# # we need to remove the periodic signal. We can do that with the savgol. You can 
# # scroll your mouse to get a a nice fit of the savgol model to the data of RS Cha

phasefold_savgol("example_data/extract_lc_data/TIC 323292655/TIC 323292655_combined_BCSAP64_65_66.txt", init_period=1.7, period_change=0.0025, savgol_window = 300)

# Now of course it helps if you know the orbital period to 5 significant diits!

phasefold_savgol("example_data/extract_lc_data/TIC 323292655/TIC 323292655_combined_BCSAP64_65_66.txt", init_period=1.66988, period_change=0.001, savgol_window = 250)

# Finally we can have a look at the pulsation light curve of RS Cha in sector 65

plot_lc("example_data/extract_lc_data/TIC 323292655/TIC 323292655_combined_BCSAP64_65_66.txtremoved_from_binned_template.txt")


# For such a quick pipeline, this did a pretty good job! We can clearly see the pulsations. It looks
# like there are some long term trends in the light curve that might be of instrumental origin.
# Ideally you can run another savgol filter to filter these out as well. Maybe a function for 
# that comes soon. In the meantime, put on your coding trousers and put it together yourself :D