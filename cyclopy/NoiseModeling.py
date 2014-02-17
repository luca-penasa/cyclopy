# -*- coding: utf-8 -*-
"""
Created on Fri Jan 31 12:00:16 2014

@author: luca
"""
from __future__ import division
import numpy as np
import scipy.optimize as optimize
from pylab import plot, interactive

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
    

def median_filter(signal, windowsize=11):
    N = len(signal)
    if (windowsize%2 == 0): #is even
        lowhwin = windowsize*0.5
        highwin = lowhwin
    else: #is odd
        lowhwin = np.floor(windowsize*0.5)
        highwin = lowhwin + 1
    
    smoothed = np.zeros(N)
    for i in np.arange(N):
        smoothed[i] = getMedian(signal, i-lowhwin, i+highwin)
        
    return smoothed

def getAR1Spectrum( alpha, ro, freqs, f_nyquist, N=None ):
    """    
    Estimate the power spectrum for an AR(1) process.
    alpha       -- is the variance of the white noise
    ro          -- is the lag-one autocorrelation coefficient    
    freqs       -- real frequencies at which to compute the spectrum [Hz]
    f_nyquist   -- nyquist freq of your data [Hz]
    N           -- number of samples in original spectrum
                   for which we are replicating. It is used for 
                   normalizing the power spectrum
       
    Returns the spectrum expected for an AR1 process with these parameters
    See e.g. Mann and Lees, 1996
    """
    #precompute these:
    if N != None:
        norm = 1.0/N
    else:
        norm = 1
    b = 1 - 2 * ro * np.cos(np.pi*freqs/f_nyquist) + ro**2
    
    return norm * alpha / b
    

def residuals(args, observed, freqs, f_nyquist, N=1, maxfreq=10):
    print (args)
    alpha, ro = args

    mask = np.ones(len(freqs), dtype=bool)
    
    if maxfreq != None:
        mask[np.argwhere(freqs > maxfreq)] = False
                
    predicted = getAR1Spectrum(alpha, ro, freqs, f_nyquist, N)
#    plot(freqs[mask], predicted[mask])
    print(np.sum((predicted[mask] - observed[mask])**2))
    return (predicted[mask] - observed[mask]) 
    

def optimizeAR1Model(obs_spectrum, freqs, f_nyquist, N, maxfreq=10, init_est = None):
    interactive(True)
#    plot(freqs, obs_spectrum)
    if init_est == None:
        init_est = (0.0, 0.5)
    
    kd,cov,infodict,mesg,ier = optimize.leastsq(residuals,
                                                init_est,
                                                args=(obs_spectrum, freqs, f_nyquist, N, maxfreq),
                                                epsfcn=0.001,           
                                                ftol=1e-16,
                                                gtol=1e-16,
                                                xtol=1e-16,
                                                maxfev=10000,
                                                full_output=True)
                                                
    return kd, mesg, ier
    
    
    
def getSpectra(series, pi=4, delta=0.01):
    import mtspec
    spectra = []
    for s in series.T:
        a, b = mtspec.mtspec(s, delta, pi)
        spectra.append(a)
        
    return spectra
    
def getCI(ci, spectra):
    import scipy.stats
    values = np.zeros(len(spectra))
    for i, s in enumerate(spectra):
        values[i] = scipy.stats.scoreatpercentile(s, ci*100)
        
    return values
        
        
        
    
