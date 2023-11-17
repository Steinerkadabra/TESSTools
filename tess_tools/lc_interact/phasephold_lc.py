import numpy as np
import matplotlib.pyplot as plt
import os
import lightkurve as lk
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter

def binning(time, mag, period, num):
    time = (time%period)/period# (1/0.5988495842998753)*0.5988495842998753
    bins = []
    means = []
    sds = []
    bins.append(0)
    ind = np.where(np.logical_and(time >= 0, time <= 0.5 / num))
    means.append(np.mean(mag[ind]))
    sds.append(np.std(mag[ind]) / np.sqrt(len(mag[ind])))
    for i in range(num - 1):
        ind = np.where(np.logical_and(time >= (i + 0.5) / num, time <= (i + 1.5) / num))
        if ind[0].size > 0:
            bins.append((i + 1) / num)
            means.append(np.mean(mag[ind]))
            sds.append(np.std(mag[ind]))
    bins.append(1)
    ind = np.where(np.logical_and(time >= (num - 0.5) / num, time <= 1))
    means.append(np.mean(mag[ind]))
    sds.append(np.std(mag[ind]))
    return np.array(bins), np.array(means), np.array(sds)



def phasefold_binning(folder, init_period = 5.0, period_change = 0.025, figsize = (10,6), invert_yaxis  = True, save_name = '', num_bins = 100,):
    count = 0
    all_files = os.listdir(folder)
    files = [file for file in all_files if file[-4:] == '.txt']
    print('Choose from the available files for this folder:')
    for file in files:
        print(count, ')', file, '\n')
        count += 1
    x = input('Which file should be plotted?')
    active_file = files[int(x)]
    data = np.loadtxt(folder + '/' + active_file).T
    lc = lk.LightCurve(data[0], data[1])



    fig, ax = plt.subplots(figsize = figsize)
    ax.plot((lc.time.value%init_period)/init_period, lc.flux.value, 'ko', ms = 1)
    ax.plot((lc.time.value%init_period)/init_period + 1, lc.flux.value, 'ko', ms = 1)
    if num_bins > 0:
        bins, means, sds = binning(lc.time.value, lc.flux.value, init_period, num_bins)
        ax.plot(bins, means, 'r' + 'o', ms=5, zorder = 5)
        ax.plot(bins + 1, means, 'r' + 'o', ms=5, zorder = 5)
    ax.set_xlabel('phase  with period {:.3f}'.format(init_period) )
    ax.set_ylabel('magnitude - mean')
    if invert_yaxis:
        ax.invert_yaxis()
    plt.tight_layout()
    def scroll_period(event):

        for line in ax.lines:
            line.remove()

        plotted_stuff.clear()
        ax.relim(visible_only=True)

        period.append(period[-1]+ period_change * event.step)

        ax.set_xlabel('phase with period {:.3f}'.format(period[-1]) )
        ax.set_ylabel('magnitude - mean')
        ax.set_xlim(0, 2)

        plotted_stuff.append(ax.plot((lc.time.value%period[-1])/period[-1], lc.flux.value, 'ko', ms = 1))
        plotted_stuff.append(ax.plot((lc.time.value%period[-1])/period[-1] + 1, lc.flux.value, 'ko', ms = 1))
        if num_bins > 0:
            bins, means, sds = binning(lc.time.value, lc.flux.value, period[-1], num_bins)
            ax.plot(bins, means, 'r' + 'o', ms=5, zorder = 5)
            ax.plot(bins + 1, means, 'r' + 'o', ms=5, zorder = 5)
        # if invert_yaxis:
        #     ax.invert_yaxis()
        plt.tight_layout()

        fig.canvas.draw()

    plotted_stuff = []
    period = [init_period]
    cid2 = fig.canvas.mpl_connect('scroll_event', scroll_period)
    plt.show()
    plt.close()
    if save_name != '':
        fig, ax = plt.subplots(figsize = figsize, dpi = 200)
        ax.plot((lc.time.value%init_period)/init_period, lc.flux.value, 'ko', ms = 1)
        ax.plot((lc.time.value%init_period)/init_period + 1, lc.flux.value, 'ko', ms = 1)
        ax.set_xlabel('phase  with period {:.3f}'.format(init_period) )
        ax.set_ylabel('magnitude - mean')
        if num_bins > 0:
            bins, means, sds = binning(lc.time.value, lc.flux.value, period[-1], num_bins)
            ax.plot(bins, means, 'r' + 'o', ms=5, zorder = 5)
            ax.plot(bins + 1, means, 'r' + 'o', ms=5, zorder = 5)
        if invert_yaxis:
            ax.invert_yaxis()
        plt.tight_layout()
        plt.savefig(folder + '/' + save_name + '.png')
        print('saved ', folder + '/' + save_name + '.png')

    fperiod  = period[-1]
    bins, means, sds = binning(lc.time.value, lc.flux.value, fperiod, num_bins)
    model = interp1d(bins, means)
    x = np.linspace(0, 1, 1000)
    fig, ax = plt.subplots()
    ax.plot(x, model(x))
    plt.show()

    ftime = lc.time.value
    fflux = lc.flux.value
    phases = ((ftime%fperiod)/fperiod)%1
    modelflux = model(phases)

    fig, ax = plt.subplots()
    ax.plot(ftime, fflux, 'ko', ms = 0.8)
    ax.plot(ftime, modelflux, 'r-', ms = 0.8)
    if invert_yaxis:
        ax.invert_yaxis()
    plt.show()

    fig, ax = plt.subplots()
    ax.plot(ftime, fflux-modelflux, 'ko', ms = 0.8)
    if invert_yaxis:
        ax.invert_yaxis()
    plt.show()
    np.savetxt(folder + '/' + active_file + 'removed_from_binned_template.txt', np.array([ftime, fflux-modelflux]).T)



def phasefold_savgol(folder, init_period = 5.0, period_change = 0.025, figsize = (10,6), invert_yaxis  = True, save_name = '', savgol_window= 300, savgol_poly= 3):
    count = 0
    all_files = os.listdir(folder)
    files = [file for file in all_files if file[-4:] == '.txt']
    print('Choose from the available files for this folder:')
    for file in files:
        print(count, ')', file, '\n')
        count += 1
    x = input('Which file should be plotted?')
    active_file = files[int(x)]
    data = np.loadtxt(folder + '/' + active_file).T
    lc = lk.LightCurve(data[0], data[1])




    fig, ax = plt.subplots(figsize = figsize)
    folded_time = (lc.time.value%init_period)/init_period
    ax.plot(folded_time, lc.flux.value, 'ko', ms = 1)
    ax.plot(folded_time + 1, lc.flux.value, 'ko', ms = 1)
    sort = np.argsort(folded_time)
    savgol = savgol_filter(lc.flux.value[sort], savgol_window, savgol_poly)
    ax.plot(folded_time[sort] + 1, savgol, 'ro', ms = 2)
    ax.plot(folded_time[sort] , savgol, 'ro', ms = 2)
    ax.set_xlabel('phase  with period {:.3f}'.format(init_period) )
    ax.set_ylabel('magnitude - mean')
    if invert_yaxis:
        ax.invert_yaxis()
    plt.tight_layout()
    def scroll_period(event):

        for line in ax.lines:
            line.remove()

        plotted_stuff.clear()
        ax.relim(visible_only=True)

        period.append(period[-1]+ period_change * event.step)

        ax.set_xlabel('phase with period {:.3f}'.format(period[-1]) )
        ax.set_ylabel('magnitude - mean')
        ax.set_xlim(0, 2)

        folded_time = (lc.time.value%period[-1])/period[-1]
        plotted_stuff.append(ax.plot(folded_time, lc.flux.value, 'ko', ms = 1))
        plotted_stuff.append(ax.plot(folded_time + 1, lc.flux.value, 'ko', ms = 1))
        

        sort = np.argsort(folded_time)
        savgol = savgol_filter(lc.flux.value[sort], savgol_window, savgol_poly)
        ax.plot(folded_time[sort] + 1, savgol, 'ro', ms = 2)
        ax.plot(folded_time[sort], savgol, 'ro', ms = 2)
        # if invert_yaxis:
        #     ax.invert_yaxis()
        plt.tight_layout()

        fig.canvas.draw()

    plotted_stuff = []
    period = [init_period]
    cid2 = fig.canvas.mpl_connect('scroll_event', scroll_period)
    plt.show()
    plt.close()
    if save_name != '':
        fig, ax = plt.subplots(figsize = figsize, dpi = 200)

        folded_time = (lc.time.value%init_period)/init_period
        ax.plot(folded_time, lc.flux.value, 'ko', ms = 1)
        ax.plot(folded_time + 1, lc.flux.value, 'ko', ms = 1)
        ax.set_xlabel('phase  with period {:.3f}'.format(init_period) )
        ax.set_ylabel('magnitude - mean')
        

        sort = np.argsort(folded_time)
        savgol = savgol_filter(lc.flux.value[sort], savgol_window, savgol_poly)
        ax.plot(folded_time[sort] + 1, savgol, 'ro', ms = 2)
        ax.plot(folded_time[sort], savgol, 'ro', ms = 2)
        if invert_yaxis:
            ax.invert_yaxis()
        plt.tight_layout()
        plt.savefig(folder + '/' + save_name + '.png')
        print('saved ', folder + '/' + save_name + '.png')

    fperiod  = period[-1]
    folded_time = (lc.time.value%fperiod)/fperiod
    sort = np.argsort(folded_time)
    savgol = savgol_filter(lc.flux.value[sort], savgol_window, savgol_poly)
    model = interp1d(folded_time[sort], savgol, bounds_error = False)
    x = np.linspace(0, 1, 1000)
    fig, ax = plt.subplots()
    ax.plot(x, model(x))
    if invert_yaxis:
        ax.invert_yaxis()
    plt.show()

    ftime = lc.time.value
    fflux = lc.flux.value
    phases = ((ftime%fperiod)/fperiod)%1
    modelflux = model(phases)

    fig, ax = plt.subplots()
    ax.plot(ftime, fflux, 'ko', ms = 0.8)
    ax.plot(ftime, modelflux, 'r-', ms = 0.8)
    if invert_yaxis:
        ax.invert_yaxis()
    plt.show()

    fig, ax = plt.subplots()
    ax.plot(ftime, fflux-modelflux, 'ko', ms = 0.8)
    if invert_yaxis:
        ax.invert_yaxis()
    plt.show()
    np.savetxt(folder + '/' + active_file + 'removed_from_binned_template.txt', np.array([ftime, fflux-modelflux]).T)