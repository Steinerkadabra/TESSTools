import numpy as np
import lightkurve as lk

def convert_to_mag(lc):
    """
    Converts a LighCurve object to magnitudes.

    :param lc: lightcurve object
    :return: reduced light curve object
    """

    lc = lc.remove_nans()

    flux = lc.flux.value
    flux_err = lc.flux_err.value #+ (np.abs(2 * np.amin(lc.flux)) if np.amin(lc.flux) < 0 else 100)
    time = lc.time.value


    flux = -2.5 * np.log10(flux/15000) + 10
    flux_err = -2.5 * np.log10(flux_err/15000) + 10

    time = time[~np.isnan(flux)]
    flux_err = flux[~np.isnan(flux)]
    flux = flux[~np.isnan(flux)]


    return lk.LightCurve(time = time, flux = flux, flux_err = flux_err)

