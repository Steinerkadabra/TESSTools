import os
import lightkurve as lk
import numpy as np
import matplotlib.pyplot as plt

from ..utils.lc_convert import convert_to_mag

def download_SC(target, sectors, combine = True, folder ='', conv_to_mag = True, remove_mean = True):
    """
    Function to download and convert SC observations.
    :param target: name of the target. Output will be saved in a file with correspodning name.
    :param sectors: array of sectors of which you want to comput the light curve.
    :param combine: if True and multiple sectors are given, combine the light curves at the end.
    :param folder: name of the output folder, default = None.
    :return:
    """
    if folder != '':
        folder = folder + '/'
    try:
        os.mkdir(folder+target)
    except FileExistsError:
        pass
    saps = []
    pdcsaps = []
    final_string = ''
    for sector in sectors:
        data = lk.search_lightcurve(target , sector = sector, author = "SPOC", cadence = "short").download()
        sap = lk.LightCurve(time = data.time, flux= data.sap_flux, flux_err= data.sap_flux_err).remove_nans()
        pdcsap = lk.LightCurve(time = data.time, flux = data.pdcsap_flux, flux_err= data.pdcsap_flux_err).remove_nans()

        if conv_to_mag:
            sap = convert_to_mag(sap).remove_outliers(4,maxiters=1)
            pdcsap = convert_to_mag(pdcsap).remove_outliers(4,maxiters=1)
        saps.append(sap)
        pdcsaps.append(pdcsap)

        if remove_mean:
            np.savetxt(folder + target + '/' + target + 'SAP' + str(sector) + '.txt',
                   np.array([sap.time.value, sap.flux.value - np.mean(sap.flux.value)]).T)
            np.savetxt(folder + target + '/' + target + 'PDCSAP' + str(sector) + '.txt',
                   np.array([pdcsap.time.value, pdcsap.flux.value - np.mean(pdcsap.flux.value)]).T)
        else:
            np.savetxt(folder + target + '/' + target + 'SAP' + str(sector) + '.txt',
                   np.array([sap.time.value, sap.flux.value]).T)
            np.savetxt(folder + target + '/' + target + 'PDCSAP' + str(sector) + '.txt',
                   np.array([pdcsap.time.value, pdcsap.flux.value]).T)

        final_string = final_string + f'{sector}_'

    if combine:
        fig, ax = plt.subplots(figsize = (10,6), dpi= 200)

        time = []
        mag = []
        for lc in saps:
            for t in range(len(lc.time)):
                time.append(lc.time.value[t])
                if remove_mean:
                    mag.append(lc.flux.value[t] - np.mean(lc.flux.value))
                else:
                    mag.append(lc.flux.value[t])

        np.savetxt(folder + target + '/' + target + 'SAP' + final_string + '.txt',
               np.array([time, mag]).T)
        ax.plot(time, mag, 'bo', ms = 0.75, label = 'SAP')
        time = []
        mag = []
        for lc in pdcsaps:
            for t in range(len(lc.time)):
                time.append(lc.time.value[t])
                if remove_mean:
                    mag.append(lc.flux.value[t] - np.mean(lc.flux.value))
                else:
                    mag.append(lc.flux.value[t])

        np.savetxt(folder + target + '/' + target + 'PDCSAP' + final_string + '.txt',
               np.array([time, mag]).T)
        ax.plot(time, mag, 'go', ms = 0.75, label = 'PDCSAP')

        ax.set_xlabel("Time")
        if convert_to_mag:
            ax.set_ylabel("Magnitude")
            ax.invert_yaxis()
        else:
            ax.set_ylabel("Flux")
        ax.legend()
        
        plt.show()



