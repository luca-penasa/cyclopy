# -*- coding: utf-8 -*-
"""
Created on Wed May 16 17:23:07 2012

@author: luca
"""

from __future__ import division

__name__ = "SignalFilters"



import numpy as np
import progressbar as pbar
import pylab
import scipy
import scipy.linalg
from matplotlib import mlab
import mtspec
from scipy import signal

pb = pbar.ProgressBar()

#Plot frequency and phase response
def mfreqz(b,a=1):
    w,h = signal.freqz(b,a) #compute freq response
    h_dB = 20 * np.log10 (abs(h)) #as log scale
    pylab.subplot(211)
    pylab.plot(w/max(w),h_dB)
    pylab.ylim(-150, 5)
    pylab.ylabel('Magnitude (db)')
    pylab.xlabel(r'Normalized Frequency (x$\pi$rad/sample)')
    pylab.title(r'Frequency response')
    pylab.subplot(212)
    h_Phase = np.unwrap(np.arctan2(np.imag(h),np.real(h))) #compute the phase
    pylab.plot(w/max(w),h_Phase)
    pylab.ylabel('Phase (radians)')
    pylab.xlabel(r'Normalized Frequency (x$\pi$rad/sample)')
    pylab.title(r'Phase response')
    pylab.subplots_adjust(hspace=0.5)

#Plot step and impulse response
def impz(b,a=1):
    l = len(b)
    impulse = np.repeat(0.,l); impulse[0] =1.
    x = np.arange(0,l)
    response = signal.lfilter(b,a,impulse)
    pylab.subplot(211)
    np.stem(x, response)
    pylab.ylabel('Amplitude')
    pylab.xlabel(r'n (samples)')
    pylab.title(r'Impulse response')
    pylab.subplot(212)
    step = np.cumsum(response)
    np.stem(x, step)
    pylab.ylabel('Amplitude')
    pylab.xlabel(r'n (samples)')
    pylab.title(r'Step response')
    pylab.subplots_adjust(hspace=0.5)

def designLinearBandpass(fa, fb, s_step, n=1001):
    f_sampling = 1/s_step    
    nyq = 0.5*f_sampling
    
    fa_ny = fa / nyq
    fb_ny = fb / nyq
    
    a = signal.firwin(n, cutoff = fa_ny, window = 'hamming')
    #Highpass filter with spectral inversion
    b = - signal.firwin(n, cutoff = fb_ny, window = 'hamming')
    b[n/2] = b[n/2] + 1
    #Combine into a bandpass filter
    d = - (a+b); d[n/2] = d[n/2] + 1
    #Frequency response
    mfreqz(d)
    pylab.show()
    
    return d

def bandpass(y, fa, fb, spacing=0.005):
	'''
	from Muller and macdonald: ice ages and astronomical causes
	'''
	f1 = np.min([fa, fb])
	f2 = np.max([fa, fb])
	n = len(y)
	
	fNyq = 1/(2*spacing)
	y = y - np.mean(y)
	ft = np.fft.fft(y)
	fre = np.linspace(0, 2*fNyq, n)
	#find band's start/end positions
	#START
	differences = np.abs(fre - f1)
	k1 = np.where(differences == np.min(differences))
	k1 = k1[0]
	if len(k1) != 1: 
		k1 = k1[0]
	#END
	differences = np.abs(fre - f2)
	k2 = np.where(differences == np.min(differences))
	k2 = k2[0]
	if len(k2) != 1: 
		k2 = k2[0]
		
	#And conjugate points
	k3 = n - 1 - k2
	k4 = n - 1 - k1
	
	#set the band's outsides to zero
	ft[0:k1] = 0
	ft[k2:k3] = 0
	ft[k4::] = 0
	
	#IFFT
	ynew = np.real(np.fft.ifft(ft))
	return ynew





def vssanan(ts, m, f, n_pc):
    '''
    Slightly-modified-python implementation of the function by:
    Schoellhamer, D. H., 2001, Singular spectrum analysis for time series with missing data:
    Geophysical Research Letters, v. 28 , no. 16 , p. 3187-3190. 
    URL: http://ca.water.usgs.gov/ja/grl/ssam.pdf
    Originally written for matlab.
    Rewritten by Luca Penasa. luca.penasa@gmail.com, may, 2012
    See http://ca.water.usgs.gov/ja/grl/ for the original one.
    
    purpose: given time series ts and window size m (m*dt), compute
             principal components xk, eigenvectors ro, and eigenvalues lam.
             see Vautard and Ghil 1989 for notation
             this version ignores nans

    inputs:  ts=time series with constant dt
             m=window size, in time steps
             f=fraction (0<f<=1) of good data points for determining pc's
             n_pc=number of reconstructed PCs to return

    outputs: xk=principal components zero mean time series, row k = pc k
             ro=eigenvectors ro(j,k), k=corresponding pc
             lam=eigenvalue vector, sorted from greatest to smallest
     
    '''
    #normalize the series
    igood = pylab.find(-np.isnan(ts)) #good indices - not-nan elements
    xmean = np.mean(ts[igood]) #mean of good indices
    xstd = np.std(ts[igood]) #std of good elements
    x = (ts - xmean) /xstd #normalized signal
    n = x.size #it's size
    
    c = np.zeros(m) #where to put correlations
    #compute autocorrelation at the various lag step 0 to m
    for j in np.arange(m):
        prod = x[0:n-j] * x[j:n]
        igood = pylab.find(-np.isnan(prod))
        c[j] = np.sum(prod[igood]) / igood.size 
        
    #now create the covariance matrix a, divide by m
    a = scipy.linalg.toeplitz(c) / m

    #compute eigenvalues / eignvectors
    lam, z = np.linalg.eig(a) #get eigenvalues and eigenvectors    
    
    #now we sort depending on lambdas from bigger to lower
    ilam = np.argsort(lam)
    ilam_flipped = np.flipud(ilam)
    lam = lam[ilam_flipped]
    
    summ = 0
    kmax = 0
    n_to_plot = np.min(np.array([m, 30])) #n of eigenvalues to plot
    pylab.semilogy(np.arange(n_to_plot) + 1,  lam[0:n_to_plot] , 'o') 
    
    ro = np.zeros((z.shape[0], n_pc))
    for i in np.arange(n_pc):
        ro[:,i] = z[:, ilam_flipped[i]]
        summ = summ + lam[i]
        print("Eigenvalue of PC " + str(i) + " is " + str(lam[i]) + ", sum of variance is " + str(summ) )
    
    
    #reconstruct the time series for the requested PCs
    xk = np.nan*np.ones(( n_pc, n - m ))
    for k in np.arange(n_pc):
        print ("Reconstructing Time series " + str(k))
        for i in np.arange(n-m):
            prod = x[i:i+m] * ro[:,k]
            igood = pylab.find(-np.isnan(prod))
            ngood = igood.size
            #this is for checking nan vs. good values, we need to respect a given ratio
            if ngood >= m * f:
                xk[k, i] = np.sum(prod[igood]) * m / ngood
    
    return xk, ro, lam


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