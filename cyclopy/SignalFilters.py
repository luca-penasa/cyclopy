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

class Spectrogram():
    def __init__(self, time_series):
        self.x_ = time_series.getPositionVector()
        self.y_ = time_series.getSeries()
        self.series_ = time_series
        self.spacing_ = time_series.getSamplingStep()
        self.detrend_method_ = 'linear'
        self.max_normalization_ = True
        self.sub_fact_ = 1
        self.highfreq_ = 10
        self.method_ = 'mtspec'
        self.win_len_ = 1
        self.N_zeros_ = 1000
        self.pi_ = 1
    
    def setWindowLenght(self, win_len = 1):
        self.win_len_ = win_len
    
    def setMethod(self, method = 'mtspec'):
        self.method_ = 'mtspec'
    
    def setHighFreq(self, highfreq=10):
        self.highfreq_ = highfreq
    
    def setSamplingStep(self, step):
        self.spacing_ = step
    
    def setZeroPadding(self, pad_to = 1000):
        self.N_zeros_ = pad_to
    
    def setSubsamplingFactor(self, sub_fact = 1):
        self.sub_fact_ = sub_fact
    
    def setMTMpi(self, pi = 1):
        self.pi_ = pi
    
    def plot(self, newfig = True, use_log=False):
        if newfig:
            pylab.figure()
        if use_log:
            pylab.imshow(np.log(self.ffts_.transpose()), interpolation='bicubic', aspect='auto', origin='Lower', extent=[min(self.eval_pos_), max(self.eval_pos_), 0, self.real_maxfreq_])
        else:
            pylab.imshow(self.ffts_.transpose(), interpolation='bicubic', aspect='auto', origin='Lower', extent=[min(self.eval_pos_), max(self.eval_pos_), 0, self.real_maxfreq_])            
            pylab.xlabel('time - ' + self.series_.unit_)
            pylab.ylabel('frequency - cycles/' + self.series_.unit_)
            pylab.colorbar()
            pylab.show()        
        
    
    
    def setLocalDetrend(self, method = 'linear'):
        '''
        method can be
        - 'linear'
        - 'mean'
        - 'none'
        '''
        self.detrend_method_ = method
        
    def setUseMaxNormalization(self, use = True):
        """
        use can be True or False
        If true each spectra is normalized placing the maximum peak at power 1
        """
        self.max_normalization_ = use
        
    def update(self):
        winlen_samples = np.int(self.win_len_/self.spacing_)
        half_winlen = np.int(winlen_samples / 2)
        print('each window have size: ' + str(half_winlen * 2))
        ffts = []
        ns = []
        
        for n in pb(np.arange(0 + half_winlen, len(self.x_) - half_winlen)[::self.sub_fact_]):
            ns.append(n)
            
            #get a piece of the signal
            small_signal = self.y_[n-half_winlen:n+half_winlen]
                        
            #perform detrending
            if self.detrend_method_ == 'linear':
                small_signal = mlab.detrend_linear(small_signal)
                
            elif self.detrend_method_ == 'mean':
                small_signal = small_signal - np.mean(small_signal)
            
            elif self.detrend_method_ == 'none': #simply go ahead
                pass
            
                           
            #do padding
            if self.N_zeros_ > len(small_signal):
                small_signal = np.concatenate((small_signal, np.zeros(self.N_zeros_-len(small_signal))))
                

            #get a spectral estimation using some method                 
            if self.method_ == 'mtspec':
                ps, f = mtspec.mtspec(small_signal, self.spacing_, self.pi_)
                
            #normalize the single spectrum
            if self.max_normalization_:
                ps = ps / np.max(ps)
    		
            #append to the ffts
            ffts.append(ps[f<=self.highfreq_])
            
            
        self.ffts_ = np.array(ffts)
        self.eval_pos_ = self.x_[ns]
        
        #the real max frequency?
        self.real_maxfreq_ = np.max(f[f<=self.highfreq_])
      
    def getSpectrum(self, position):
        '''
        get the spectrum that were computeted for the position nearest to position
        '''
        nearest = pylab.find( np.min(np.abs(self.eval_pos_ - position))  == np.abs(self.eval_pos_ - position))
        distance = np.abs(self.eval_pos_[nearest[0]] - position)
        print("Real positio for requested spectrum is: " + str(self.eval_pos_[nearest[0]]) + " with a distance of: "+ str(distance))
        return self.ffts_[nearest[0]] 
       
        

def spectrogram_mtspec(pos, signal, win_len=1.2, pi=2, spacing=0.005, highfreq=8, pad_to=1000, sub_fact=10, plot=True, max_normalize=True, use_log=False):
    winlen_samples = np.int(win_len/spacing) #number of samples into a sample signal
    half_winlen = np.int(winlen_samples / 2)
    print('each window have size: ' + str(half_winlen * 2))
    ffts = []
    ns = []
    for n in pb(np.arange(0 + half_winlen, len(signal) - half_winlen)[::sub_fact]):
        ns.append(n)
        #working on signal
        small_signal = signal[n-half_winlen:n+half_winlen]
        from matplotlib import mlab
        #small_signal = mlab.detrend_linear(small_signal)
        #small_signal = small_signal * np.hanning(len(small_signal))
        #small_signal = np.concatenate((small_signal, np.zeros(pad_to-len(small_signal))))
        #getFFT
        import mtspec
        ps, f = mtspec.mtspec(small_signal, spacing, pi, pad_to)
        if max_normalize:
            ps = ps / np.max(ps)
		
        ffts.append(ps[f<=highfreq])
        
        
    ffts = np.array(ffts)
    eval_pos = pos[ns]
    real_maxfreq = np.max(f[f<=highfreq])
    if plot == True:
        if use_log == True:
            pylab.imshow(np.log(ffts.transpose()), interpolation='bicubic', aspect='auto', origin='Lower', extent=[min(eval_pos), max(eval_pos), 0, real_maxfreq])
        else:
            pylab.imshow(ffts.transpose(), interpolation='bicubic', aspect='auto', origin='Lower', extent=[min(eval_pos), max(eval_pos), 0, real_maxfreq])            
        pylab.xlabel('position')
        pylab.ylabel('frequency - cycles/m')
        pylab.colorbar()
        pylab.show()
          
    return eval_pos, ffts
    

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
        
