# -*- coding: utf-8 -*-
"""
Created on Mon Jun 18 16:02:02 2012

@author: luca
"""

__name__ = "Spectrogram"

import numpy as np
import progressbar as pbar
import pylab

from matplotlib import mlab
import mtspec



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
        pb = pbar.ProgressBar()
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
    
