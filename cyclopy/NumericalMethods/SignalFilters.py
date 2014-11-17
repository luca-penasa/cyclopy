# -*- coding: utf-8 -*-
"""
Created on Wed May 16 17:23:07 2012

@author: luca
"""

from __future__ import division
from numpy.lib.type_check import nan_to_num

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

def designLinearBandpass(fa, fb, s_step, n=1001, show=False):
    f_sampling = 1/s_step    
    nyq = f_sampling / 2
    
    fa_ny = fa / nyq
    fb_ny = fb / nyq
    
    a = signal.firwin(n, cutoff = fa_ny, window = 'blackmanharris')
    #Highpass filter with spectral inversion
    b = - signal.firwin(n, cutoff = fb_ny, window = 'blackmanharris')
    # b[n/2] = b[n/2] + 1
    #Combine into a bandpass filter
    # d = - (a+b); d[n/2] = d[n/2] + 1


    b[nyq/2] = b[nyq/2] + 1
    #Combine into a bandpass filter
    d = - (a+b)
    d[nyq/2] = d[nyq/2] + 1


    #Frequency response
    if show:
        mfreqz(d)
        pylab.show()
    
    return d

def bandpass(y, fa, fb, spacing=1):
    '''
    from Muller and macdonald: ice ages and astronomical causes
    '''
    f1 = np.min([fa, fb])
    f2 = np.max([fa, fb])
    n = len(y)

    fNyq = 1/(2*spacing)
    # import scipy.signal
    # y = scipy.signal.detrend(y)

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
    k3 = n + 1 - k2
    k4 = n + 1 - k1

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

def fill_nan_linear(y, x=None):
    if x == None:
        x = np.arange(len(y))

    good_ids = np.isfinite(y)
    nan_ids = np.isnan(y)

    real_x = x[good_ids]
    real_y = y[good_ids]

    from scipy.interpolate import interp1d
    interpolator  = interp1d(real_x, real_y, 'linear')

    new_y = y[:]
    new_y[nan_ids] = interpolator(x[nan_ids])

    return new_y

def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    r"""Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techniques.
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.
    Examples
    --------
    t = np.linspace(-4, 4, 500)
    y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    ysg = savitzky_golay(y, window_size=31, order=4)
    import matplotlib.pyplot as plt
    plt.plot(t, y, label='Noisy signal')
    plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    plt.plot(t, ysg, 'r', label='Filtered signal')
    plt.legend()
    plt.show()
    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """
    import numpy as np
    from math import factorial

    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError, msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')

def median_filter(signal, windowsize=11):
    from cyclopy.Utils import getMedian
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