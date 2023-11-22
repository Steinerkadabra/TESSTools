import lightkurve as lk 
import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from matplotlib.colorbar import Colorbar
from matplotlib import patches
import astropy.visualization as stretching
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.io import ascii
import os
from astropy.table import Table


from ..utils.lc_utils import convert_to_mag
from ..utils.lc_utils import cut_lc
from ..utils.gaia_utils import get_gaia_data_from_tic
from ..utils.gaia_utils import add_gaia_figure_elements
from ..utils.gaia_utils import plot_orientation

def get_lc_from_FFI(target, tic, sectors, maglim = 6, cut_lcs = None,do_correction = True, combine_sectors = True, folder = '', plot_final = True, remove_mean = True, conv_to_mag = True):
    """
    Interactive plot to choose the aperture for which to calculate a light curve. Output will be the raw light curve, a Principal component analysis
    corrected lightcure and/or a simple background  corrected light curve for each sector and if activated the combined light curve of all sectors.
    Also saves the plot of the FFI and the chosen apertures.
    A lot of the code is taken from tpfplotter, especially all the Gaia stuff: https://github.com/jlillo/tpfplotter.
    :param target: name of the target. Output will be saved in a file with correspodning name.
    :param tic: TIC ID of the target.
    :param sectors: array of sectors of which you want to comput the light curve.
    :param maglim: magnitude limit for which stars to include in the plot.
    :param cut_lcs: parts of the light curve that should be cut out.
    :param do_correction: if True, apply prinicipal component analysis corrections.
    :param combine_sectors: if True and multiple sectors are given, combine the light curves at the end.
    :param folder:  name of the output folder, default = None.
    :param plot_final: if True and combine_sectors True: plot the final combined corrected light curves.
    :return:
    """

    if folder != '':
        folder = folder + '/'
    gaia_id, mag = get_gaia_data_from_tic(tic)

    BCSAP_s = []
    corrected_s = []
    raw_s = []


    tpfs = []
    print('Downloading data:')

    for sector in sectors:
        print('Sector ' , sector)
        # tpf = lk.search_targetpixelfile(target, sector=sector, exptime= "short")
        # tpf = tpf.download(quality_bitmask='default')
        # try:
        #     print("SC Cutout: ", tpf.shape[1],tpf.shape[2])
        #     tpfs.append(tpf)
        # except:
            # print('No Shortcadence D/ata available, get FFI. ')
        # print(sector)
        tpf = lk.search_tesscut(target, sector=sector)
        tpf = tpf.download(quality_bitmask='default', cutout_size=(20, 20))
        tpfs.append(tpf)



    print('Finished Downloading. Start extracting aperture.')
    sector_count = 0

    for tpf_gaia in tpfs:
        nans = []
        # print(tpf_gaia.flux_err)
        print(np.argwhere(np.isnan(tpf_gaia.flux_err)))
        for shit in np.argwhere(np.isnan(tpf_gaia.flux_err)):
            print(shit, shit[0])
            if shit[0] not in nans:
                nans.append(shit[0])
        print(np.argwhere(tpf_gaia.flux_err <= 0))
        for shit in np.argwhere(tpf_gaia.flux_err <= 0):
            print(shit, shit[0])
            if shit[0] not in nans:
                nans.append(shit[0])

        print("remove cadences:", nans)
        take = np.ones(len(tpf_gaia.flux_err), dtype=bool)
        for shit in nans:
            take[shit] = False



        tpf_gaia = tpf_gaia[take]
        aperture = np.zeros((tpf_gaia.shape[1],tpf_gaia.shape[2]), dtype=bool)
        simple_backgorund_aper = None
        simple_backgorund_aper = np.zeros((tpf_gaia.shape[1],tpf_gaia.shape[2]), dtype=bool)
        print(sectors[sector_count])
        do_BCSAP = True
        gaia_fig = plt.figure(figsize=(6.93, 5.5))
        gs = gridspec.GridSpec(1,3, height_ratios=[1], width_ratios=[1,0.05,0.01])
        gs.update(left=0.05, right=0.95, bottom=0.12, top=0.95, wspace=0.01, hspace=0.03)
        ax1 = plt.subplot(gs[0,0])
        mean_tpf = np.mean(tpf_gaia.flux.value, axis=0)
        nx, ny = np.shape(mean_tpf)
        norm = ImageNormalize(stretch=stretching.LogStretch())
        division = int(np.log10(np.nanmax(tpf_gaia.flux.value)))
        splot = plt.imshow(np.nanmean(tpf_gaia.flux.value, axis=0) / 10 ** division, norm=norm, extent=[tpf_gaia.column, tpf_gaia.column + ny, tpf_gaia.row, tpf_gaia.row + nx], origin='lower', zorder=0)


            # Gaia sources
        r, res = add_gaia_figure_elements(tpf_gaia, magnitude_limit=mag + float(maglim), targ_mag=mag)
        x, y, gaiamags = r
        x, y, gaiamags = np.array(x) + 0.5, np.array(y) + 0.5, np.array(gaiamags)
        size = 128.0 / 2 ** ((gaiamags - mag))
        plt.scatter(x, y, s=size, c='red', alpha=0.6, edgecolor=None, zorder=10)

        # Gaia source for the target
        this = np.where(np.array(res['Source']) == int(gaia_id))[0]
        plt.scatter(x[this], y[this], marker='x', c='white', s=32, zorder=11)

        # Legend
        add = 0
        if int(maglim) % 2 != 0:
            add = 1
        maxmag = int(maglim) + add
        legend_mags = np.linspace(-2, maxmag, int((maxmag + 2) / 2 + 1))
        fake_sizes = mag + legend_mags  # np.array([mag-2,mag,mag+2,mag+5, mag+8])
        for f in fake_sizes:
            size = 128.0 / 2 ** ((f - mag))
            plt.scatter(0, 0, s=size, c='red', alpha=0.6, edgecolor=None, zorder=10,
                        label=r'$\Delta m=$ ' + str(int(f - mag)))

        ax1.legend(fancybox=True, framealpha=0.7, fontsize ='x-small')

        # Source labels
        # print(this)
        # print(x, y)
        dist = np.sqrt((x - x[this]) ** 2 + (y - y[this]) ** 2)
        dsort = np.argsort(dist)
        for d, elem in enumerate(dsort):
            if dist[elem] < 6:
                plt.text(x[elem] + 0.1, y[elem] + 0.1, str(d + 1), color='white', zorder=100)

            # Orientation arrows
        plot_orientation(tpf_gaia)
        # Labels and titles
        plt.xlim(tpf_gaia.column, tpf_gaia.column + ny)
        plt.ylim(tpf_gaia.row, tpf_gaia.row + nx)
        plt.xlabel('Pixel Column Number', fontsize=16)
        plt.ylabel('Pixel Row Number', fontsize=16)

        # Colorbar
        cbax = plt.subplot(gs[0, 1])  # Place it where it should be.
        pos1 = cbax.get_position()  # get the original position
        pos2 = [pos1.x0 - 0.05, pos1.y0, pos1.width, pos1.height]
        cbax.set_position(pos2)  # set a new position

        cb = Colorbar(ax=cbax, mappable=splot, orientation='vertical', ticklocation='right')
        plt.xticks(fontsize=14)
        exponent = r'$\times 10^' + str(division) + '$'
        cb.set_label(r'Flux ' + exponent + r' (e$^-$)', labelpad=10, fontsize=16)

        dist = np.sqrt((x - x[this]) ** 2 + (y - y[this]) ** 2)
        GaiaID = np.array(res['Source'])
        srt = np.argsort(dist)
        x, y, gaiamags, dist, GaiaID = x[srt], y[srt], gaiamags[srt], dist[srt], GaiaID[srt]

        aperture_patches = []
        for i in range(aperture.shape[0]):
            for j in range(aperture.shape[1]):
                if aperture[i, j]:
                    aperture_patches.append(ax1.add_patch(patches.Rectangle((j + tpf_gaia.column, i + tpf_gaia.row),
                                                    1, 1, color='tomato', fill=True, alpha=0.4)))
                    aperture_patches.append(ax1.add_patch(patches.Rectangle((j + tpf_gaia.column, i + tpf_gaia.row),
                                                    1, 1, color='tomato', fill=False, alpha=1, lw=2)))
        if do_BCSAP:
            aperture_patches_BCSAP = []
            for i in range(simple_backgorund_aper.shape[0]):
                for j in range(simple_backgorund_aper.shape[1]):
                    if simple_backgorund_aper[i, j]:
                        aperture_patches_BCSAP.append(
                            ax1.add_patch(patches.Rectangle((j + tpf_gaia.column, i + tpf_gaia.row),
                                                            1, 1, color='magenta', fill=True,
                                                            alpha=0.4)))
                        aperture_patches_BCSAP.append(
                            ax1.add_patch(patches.Rectangle((j + tpf_gaia.column, i + tpf_gaia.row),
                                                            1, 1, color='magenta', fill=False, alpha=1,
                                                            lw=2)))

        def onclick_aperture(event):
            # print('button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
            #       (event.button, event.x, event.y, event.xdata, event.ydata))

            k = int(event.xdata) - tpf_gaia.column
            l = int(event.ydata) - tpf_gaia.row
            if event.button == 1:
                for p in aperture_patches:
                    p.remove()
                aperture_patches.clear()

                aperture[l,k] = not aperture[l,k]
                for i in range(aperture.shape[0]):
                    for j in range(aperture.shape[1]):
                        if aperture[i, j]:
                            aperture_patches.append(ax1.add_patch(patches.Rectangle((j + tpf_gaia.column, i + tpf_gaia.row),
                                                                                    1, 1, color='tomato', fill=True,
                                                                                    alpha=0.4)))
                            aperture_patches.append(ax1.add_patch(patches.Rectangle((j + tpf_gaia.column, i + tpf_gaia.row),
                                                                                    1, 1, color='tomato', fill=False, alpha=1,
                                                                                    lw=2)))

            elif event.button == 3 and do_BCSAP:
                for p in aperture_patches_BCSAP:
                    p.remove()
                aperture_patches_BCSAP.clear()
                for i in range(simple_backgorund_aper.shape[0]):
                    for j in range(simple_backgorund_aper.shape[1]):
                        if i == l and j == k and not simple_backgorund_aper[i, j] :
                            simple_backgorund_aper[i, j] = True
                        else:
                            simple_backgorund_aper[i, j] = False

                        if simple_backgorund_aper[i, j]:
                            aperture_patches_BCSAP.append(ax1.add_patch(patches.Rectangle((j + tpf_gaia.column, i + tpf_gaia.row),
                                                                                    1, 1, color='magenta', fill=True,
                                                                                    alpha=0.4)))
                            aperture_patches_BCSAP.append(ax1.add_patch(patches.Rectangle((j + tpf_gaia.column, i + tpf_gaia.row),
                                                                                    1, 1, color='magenta', fill=False, alpha=1,
                                                                                    lw=2)))

            gaia_fig.canvas.draw()

        if do_BCSAP:
            ax1.set_title('Left click: Choose aperture         Right click: choose pixel for BCSAP')
        else:
            ax1.set_title('Left click: Choose aperture')


        cid = gaia_fig.canvas.mpl_connect('button_press_event', onclick_aperture)
        plt.show()

        if True not in simple_backgorund_aper:
            do_BCSAP = False
        ax1.set_title('')

        IDs = np.arange(len(x)) + 1
        inside = np.zeros(len(x))

        for i in range(aperture.shape[0]):
            for j in range(aperture.shape[1]):
                if aperture[i, j]:
                    xtpf, ytpf = j + tpf_gaia.column, i + tpf_gaia.row
                    _inside = np.where((x > xtpf) & (x < xtpf + 1) &
                                       (y > ytpf) & (y < ytpf + 1))[0]
                    inside[_inside] = 1

        try:
            os.mkdir(folder+target)
        except FileExistsError:
            pass

        data = Table([IDs, GaiaID, x, y, dist, dist * 21., gaiamags, inside.astype('int')], names=['# ID', 'GaiaID', 'x', 'y', 'Dist_pix', 'Dist_arcsec', 'Gmag', 'InAper'])
        ascii.write(data, folder+target+ '/'+target +'stars_around_S' + str(sectors[sector_count]) + '.txt', overwrite=True)
        gaia_fig.savefig(folder+target+ '/'+target +'stars_around_S' + str(sectors[sector_count]) + '.pdf')
        plt.close()



        # for tpf in tpfs:
        tpf =tpf_gaia
        raw = tpf.to_lightcurve(aperture_mask=aperture)
        if do_correction:

            regressors = tpf.flux.value[:, ~aperture]
            dm = lk.DesignMatrix(regressors, name='regressors')
            dm = dm.pca(5)
            dm = dm.append_constant()

            corrector = lk.RegressionCorrector(raw)
            corrected_lc = corrector.correct(dm)
            corrector.diagnose()
            plt.show()
            plt.close()

            model = corrector.model_lc
            model -= np.percentile(model.flux, 5)
            corrected_lc = raw - model



        if do_BCSAP:
            raw_back = tpf.to_lightcurve(aperture_mask=simple_backgorund_aper)
            BCSAP = lk.LightCurve(raw.time.value, raw.flux.value - np.sum(aperture) * raw_back.flux.value)
            if isinstance(cut_lcs, np.ndarray):
                for cut in cut_lcs:
                    BCSAP = cut_lc(BCSAP, cut)
            if conv_to_mag:
                BCSAP = convert_to_mag(BCSAP).remove_outliers(4,maxiters=1)
        if isinstance(cut_lcs, np.ndarray):
            for cut in cut_lcs:
                print(cut)
                raw = cut_lc(raw, cut)
                if do_correction:
                    corrected_lc = cut_lc(corrected_lc, cut)
        if conv_to_mag:
            raw = convert_to_mag(raw).remove_outliers(4,maxiters=1)
        if do_correction:
            if conv_to_mag:
                corrected_lc = convert_to_mag(corrected_lc).remove_outliers(4,maxiters=1)


        fig, ax = plt.subplots()
        ax.plot(raw.time.value, raw.flux.value, 'ro', ms =1, label = 'raw')
        if do_correction:
            ax.plot(corrected_lc.time.value, corrected_lc.flux.value, 'go', ms =1, label = 'corrected')

        if do_BCSAP:
            ax.plot(BCSAP.time.value, BCSAP.flux.value, 'bo', ms =1, label = 'BCSAP')
        ax.set_xlabel('Time')
        ax.set_ylabel('Magnitude')
        ax.legend()
        if conv_to_mag:
            ax.invert_yaxis()

        plt.show()

        if do_BCSAP:
            BCSAP_s.append(BCSAP)
        if do_correction:
            corrected_s.append(corrected_lc)
        raw_s.append(raw)
        plt.close()
        sector_count += 1


    if len(sectors) == 1:
        if remove_mean:
            if do_BCSAP:
                np.savetxt(folder+target+ '/'+target +'BCSAP_' + str(sectors[0]) + '.txt', np.array([BCSAP_s[0].time.value,BCSAP_s[0].flux.value-np.mean(BCSAP_s[0].flux.value) ]).T)
            np.savetxt(folder+target+ '/'+target +'raw' + str(sectors[0]) + '.txt', np.array([raw_s[0].time.value,raw_s[0].flux.value-np.mean(raw_s[0].flux.value) ]).T)
            if do_correction:
                np.savetxt(folder+target+ '/'+target +'corrected' + str(sectors[0]) + '.txt', np.array([corrected_s[0].time.value,corrected_s[0].flux.value-np.mean(corrected_s[0].flux.value) ]).T)
        else:
            if do_BCSAP:
                np.savetxt(folder+target+ '/'+target +'BCSAP_' + str(sectors[0]) + '.txt', np.array([BCSAP_s[0].time.value,BCSAP_s[0].flux.value ]).T)
            np.savetxt(folder+target+ '/'+target +'raw' + str(sectors[0]) + '.txt', np.array([raw_s[0].time.value,raw_s[0].flux.value]).T)
            if do_correction:
                np.savetxt(folder+target+ '/'+target +'corrected' + str(sectors[0]) + '.txt', np.array([corrected_s[0].time.value,corrected_s[0].flux.value ]).T)

        if plot_final:
            time = corrected_s[0].time.value
            mag = corrected_s[0].flux.value
            if do_BCSAP:
                time_bcsap = BCSAP_s[0].time.value
                mag_bcsap = BCSAP_s[0].flux.value
            ffig, fax = plt.subplots(3, 1, figsize=(10, 8))
            fax[0].plot(time, mag, 'ko', ms=0.75)
            fax[0].set_xlabel('Time - 2457000 [BTJD days]')
            fax[0].set_ylabel('magnitude - mean')

            if do_BCSAP:
                pdg_bcsap = lk.LightCurve(time_bcsap, mag_bcsap).to_periodogram()
                fax[0].plot(time_bcsap, mag_bcsap, 'ro', ms=0.75)
                fax[2].plot(pdg_bcsap.frequency, (1000 * pdg_bcsap.power) ** 2 / pdg_bcsap.frequency, 'r-', lw=1)
                fax[1].plot(pdg_bcsap.frequency, 1000 * pdg_bcsap.power, 'r-', lw=1)
            pdg = lk.LightCurve(time, mag).to_periodogram()
            fax[1].plot(pdg.frequency, 1000 * pdg.power, 'k-', lw=1)
            fax[1].set_xlabel('frequency $d^{-1}$')
            fax[1].set_ylabel('amplitude (mmag)')
            fax[2].plot(pdg.frequency, (1000 * pdg.power) ** 2 / pdg.frequency, 'k-', lw=1)
            fax[2].set_xlabel('frequency $d^{-1}$')
            fax[2].set_ylabel('PSD')
            fax[2].set_yscale('log')
            fax[2].set_xscale('log')

            if conv_to_mag:
                fax[0].invert_yaxis()
            plt.show()

    else:
        count = 0
        for sector in sectors:
            if remove_mean:
                if do_BCSAP:
                    np.savetxt(folder+target+ '/'+target +'_BCSAP_' + str(sectors[count]) + '.txt', np.array([BCSAP_s[count].time.value,BCSAP_s[count].flux.value-np.mean(BCSAP_s[count].flux.value) ]).T)
                np.savetxt(folder+target+ '/'+target +'_raw' + str(sectors[count]) + '.txt', np.array([raw_s[count].time.value,raw_s[count].flux.value-np.mean(raw_s[count].flux.value) ]).T)
                if do_correction:
                    np.savetxt(folder+target+ '/'+target +'_corrected' + str(sectors[count]) + '.txt', np.array([corrected_s[count].time.value,corrected_s[count].flux.value-np.mean(corrected_s[count].flux.value) ]).T)
                count += 1
            else:
                if do_BCSAP:
                    np.savetxt(folder+target+ '/'+target +'_BCSAP_' + str(sectors[count]) + '.txt', np.array([BCSAP_s[count].time.value,BCSAP_s[count].flux.value  + np.mean(BCSAP_s[0].flux.value) - np.mean(BCSAP_s[count].flux.value) ]).T)
                np.savetxt(folder+target+ '/'+target +'_raw' + str(sectors[count]) + '.txt', np.array([raw_s[count].time.value,raw_s[count].flux.value + np.mean(raw_s[0].flux.value) -  np.mean(raw_s[count].flux.value) ]).T)
                if do_correction:
                    np.savetxt(folder+target+ '/'+target +'_corrected' + str(sectors[count]) + '.txt', np.array([corrected_s[count].time.value,corrected_s[count].flux.value+ np.mean(corrected_s[0].flux.value) -  np.mean(corrected_s[count].flux.value) ]).T)
                count += 1

        if combine_sectors:
            time = []
            mag = []
            mean_0 = np.mean(corrected_s[0].flux.value)
            for lc in corrected_s:
                for t in range(len(lc.time.value)):
                    time.append(lc.time.value[t])
                    if remove_mean:
                        mag.append(lc.flux.value[t] - np.mean(lc.flux.value))
                    else:
                        mag.append(lc.flux.value[t]- np.mean(lc.flux.value) + mean_0)

            string = str(sectors[0])
            for t in range(1, len(sectors)):
                string = string + '_' + str(sectors[t])

            np.savetxt(folder + target + '/' + target + '_combined_corrected' + string + '.txt',
                   np.array([time, mag]).T)
            if plot_final:
                ffig , fax = plt.subplots(3,1, figsize = (10,8))
                fax[0].plot(time, mag, 'ko', ms = 0.75)
                fax[0].set_xlabel('Time - 2457000 [BTJD days]')
                fax[0].set_ylabel('magnitude - mean')
                pdg = lk.LightCurve(time, mag).to_periodogram()
                fax[1].plot(pdg.frequency, 1000*pdg.power, 'k-', lw = 1)
                fax[1].set_xlabel('frequency $d^{-1}$')
                fax[1].set_ylabel('amplitude (mmag)')
                fax[2].plot(pdg.frequency, (1000*pdg.power)**2/pdg.frequency, 'k-', lw = 1)
                fax[2].set_xlabel('frequency $d^{-1}$')
                fax[2].set_ylabel('PSD')
                fax[2].set_yscale('log')
                fax[2].set_xscale('log')

                if conv_to_mag:
                    fax[0].invert_yaxis()
                plt.show()
                plt.close()

            if do_BCSAP:
                time = []
                mag = []
                mean_0 = np.mean(BCSAP_s[0].flux.value)
                for lc in BCSAP_s:
                    for t in range(len(lc.time.value)):
                        time.append(lc.time.value[t])
                        if remove_mean:
                            mag.append(lc.flux.value[t]- np.mean(lc.flux.value))
                        else:
                            mag.append(lc.flux.value[t]- np.mean(lc.flux.value) + mean_0)

                string = str(sectors[0])
                for t in range(1, len(sectors)):
                    string = string + '_' + str(sectors[t])

                np.savetxt(folder + target + '/' + target + '_combined_BCSAP' + string + '.txt',
                       np.array([time, mag]).T)