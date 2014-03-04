# -*- coding: utf-8 -*-
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

import progressbar as pbar

pb = pbar.ProgressBar()

__name__ = 'SpectralEstimation'
__author__ = 'Luca Penasa'


def mtm_from_R(x, K=3, NW=6, nFFT='default', plot=0, deltat=0.001, jkCIProb = 0.95, maxf=10, plotftest=True):
    from rpy2 import robjects
    from rpy2.robjects.packages import importr
    from rpy2.robjects import r
    MTM = importr('multitaper')
    #creating the R-type serie
    serie = robjects.FloatVector(x)
    serie = r.ts(serie, deltat = deltat)

    #executing the mtm function
    results = MTM.spec_mtm(serie,k=K, nw=NW, nFFT=nFFT, plot=0, Ftest=1, jackknife=1, jkCIProb = jkCIProb)
    #extracting variables from results
    freqs = np.array(results.rx2('freq'))
    spec = np.array(results.rx2('spec'))
    upperCI = np.array(results.rx2('mtm').rx2('jk').rx2('upperCI'))
    lowerCI = np.array(results.rx2('mtm').rx2('jk').rx2('lowerCI'))
    Ftest = np.array(results.rx2('mtm').rx2('Ftest'))
    if plot == 1:
        fig = plt.figure()
        plt.hold(1)
        ax1 = fig.add_subplot(111)
        ax1.plot(freqs, spec, 'k-', linewidth=2)
        ax1.plot(freqs, upperCI, 'r:')
        ax1.plot(freqs, lowerCI, 'g:')
        ax1.set_ylabel('Power Spectral Density')
        ax1.set_xlabel('Frequency [cycles/m]')

        ax2 = ax1.twinx()
        if plotftest == True:
            ax2.plot(freqs, Ftest, 'y')
            ax2.set_ylabel('Ftest')
            ax2.set_xlim(0, maxf)
        else:
            ax1.set_xlim(0, maxf)
        plt.show()

    return freqs, spec, upperCI, lowerCI, Ftest



