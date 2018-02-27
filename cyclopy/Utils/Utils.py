from __future__ import division

from cyclopy.External import peakdetect
from copy import deepcopy

import numpy as np


from matplotlib.pyplot import plot, text, title, fill_between, ylabel, ylim, xlim, grid, xlabel, ylabel, gca
from cyclopy.NumericalMethods import smooth_signal


def getMedian(signal, startarg, stoparg):
    """
    get the median value of a slice of an array
    stoparg is excluded
    startarg comprised
    """
    if startarg < 0:
        startarg = 0

    if stoparg > len(signal):
        stoparg = len(signal)

    return np.median(signal[startarg:stoparg])


def get_numbered_filename(filename, extension, path="./"):
    """
    get a suitable filename for saving a file.
    it will return a complete filename of type:
    [filename]_0.[extension] with numbering chosen so to not overwrite existing files in [path]

    the path wil be appended tot the returned filename so that the retuning string may be used to save the file directly
    """
    import os.path

    i = 0

    while os.path.isfile(path + "/" + filename + "_" + str(i) + "." + extension):
        i += 1

    filename = path + "/" + filename + "_" + str(i) + "." + extension


    return filename


def annotate_peaks(f, powers, noise=None, dots=False, fontsize=8, lookahead=10, **kwargs):
    """
    annotate_peaks plot periods of peaks in current figure.
    noise is the expected noise level for each spctra power and it is used to extract peaks

    all args should be numpy 1d arrays
    """

    newp = deepcopy(powers)

    if noise is not None:
        newp[powers <= noise] = 0  # put to zero powers lowers than noise


    # now we need to locate peaks

    max_ids, min_ids = peakdetect(powers, lookahead=lookahead)  # we detect peaks on the original powers

    if max_ids is None:
        raise IndexError("no peaks detected!")

    max_ids = [i[0] for i in max_ids if newp[i[0]] != 0]
    # min_ids = [i[0] for i in min_ids if newp[i[0]] != 0]



    # print(max_ids)
    if dots:
        plot(f[max_ids], powers[max_ids], 'o')

    # put labels over the peaks


    maxp = np.max(powers)
    overhang = 0.01 * maxp

    for id in max_ids:
        text(f[id], powers[id] + overhang, str(round(1 / f[id], 2)), rotation=90, horizontalalignment='center',
             verticalalignment='bottom', fontsize=fontsize, **kwargs)

    return max_ids


def plot_spectrum_with_noise(freq, spectrum, freq2=None, noise_spectrum=None, plottitle="spectrum",
                             smoothed=None, infigure_text=None, maxy=None, gridon=False,
                             xlabeltext="Frequency cycles / time", ylabeltext="Power",
                             spectrum_color=[1, 0, 0, 0.8], noise_color=[0, 1, 0, 0.4],
                             fontsize=8):
    """
    notice it will plot in the currecnt active figure/subplot
    smoothed        : overlay a smoothed version of the spectrum with computed via gaussian smoothing with bandwidth smoothed

    """
    title(plottitle)  # title of this plot

    #create a smoothed-out version of the signal if requested
    if smoothed is not None:
        smoothed = smooth_signal(freq, spectrum, smoothed)
        plot(freq, smoothed)

    if infigure_text is not None:
        x = np.max(freq) - 0.1 * np.max(freq)
        y = np.max(spectrum) - 0.1 * np.max(spectrum)

        text(x, y, infigure_text, fontdict={'fontsize': fontsize})



    # plot the power
    fill_between(freq, spectrum, facecolor=spectrum_color, edgecolor=spectrum_color)

    # add the noise if given
    if (freq2 is not None) & (noise_spectrum is not None):
        fill_between(freq2, noise_spectrum, facecolor=noise_color, edgecolor=noise_color)

    if maxy is None:
        if noise_spectrum is not None:
            maxy = np.max([np.max(spectrum), np.max(noise_spectrum)])
        else:
            maxy = np.max(spectrum)

        maxy += 0.1 * maxy

    ylim(0, maxy)

    xlim(0, np.max(freq))

    if gridon:
        grid()

    annotate_peaks(freq, spectrum, noise_spectrum, fontsize=fontsize)
    # legend()

    xlabel(xlabeltext)
    ylabel(ylabeltext)

    ax = gca()
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    ax.tick_params(axis='both', direction='out')
    ax.get_xaxis().tick_bottom()  # remove unneeded ticks
    ax.get_yaxis().tick_left()


def generateSinTS(period = 100, delta=1, length=1000): #e.g kyrs
    N = np.floor(length / delta)
    from ..Elements import TimeSeriesEven
    x = np.arange(N)  / delta
    y = np.sin(x/period * 2*np.pi)
    return TimeSeriesEven(y, delta, 0)

