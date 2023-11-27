import os
import matplotlib.pyplot as plt
import lightkurve as lk
import numpy as np


def plot_lc(file, figsize = (10,6), invert_yaxis  = True):
    count = 0
    if os.path.isfile(file):
        data = np.loadtxt(file).T
    else:
        folder = file
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
    pdg = lc.to_periodogram()
    fig, ax = plt.subplots(2, 1, figsize = figsize)
    ax[0].plot(lc.time.value, lc.flux.value, 'ko', ms = 1)
    if invert_yaxis:
        ax[0].invert_yaxis()
        ax[0].set_ylabel('magnitude - mean')
    else:
        ax[0].set_ylabel('flux - mean')
    ax[1].plot(pdg.frequency.value, pdg.power.value, 'k-')
    ax[0].set_xlabel('time in days')
    ax[1].set_xlabel('frequency in 1/d')
    ax[1].set_ylabel('power')
    plt.tight_layout()
    plt.show()
    plt.close()