# -*- coding: utf-8 -*-
from __future__ import division 
import numpy as np
from matplotlib import mlab
import pylab 
import matplotlib.pyplot as plt

import progressbar as pbar

pb = pbar.ProgressBar()

__name__ = 'SpectralEstimation'
__author__ = 'Luca Penasa'

#luca.penasa[at]gmail.com

def spectrogram(pos, signal, win_len=1.2, spacing=0.001, highfreq=15, pad_to=5000, plot=True):
	winlen_samples = int(win_len/spacing)
	maxfreqnumber = int(highfreq * spacing * pad_to)
	ffts = []
	ns = []
	for n in np.arange(0, len(signal) - winlen_samples):
		ns.append(n)
		#working on signal
		small_signal = signal[n:n+winlen_samples]
		small_signal = mlab.detrend_linear(small_signal)
		small_signal = small_signal * np.hanning(len(small_signal))
		small_signal = np.concatenate((small_signal, np.zeros(pad_to-len(small_signal))))
		#getFFT
		ps = np.abs(np.fft.fft(small_signal)[0:maxfreqnumber])**2
		ffts.append(ps / max(ps))
	ffts = np.array(ffts)
	eval_pos = pos[ns]

	if plot == True:
		pylab.imshow(ffts.transpose(), interpolation='bilinear', aspect='auto', origin='Lower', extent=[min(eval_pos), max(eval_pos), 0, highfreq])
		pylab.xlabel('position - cm')
		pylab.ylabel('frequency - cycles/m')
		pylab.colorbar()
		pylab.show()
	return ffts

def spectrogram_ok(x, signal, win_len=4, maxfreq = 8, pad_to=1000, samples_spacing = 0.1, pi = 2, plot=True ):
	import mtspec
	#x must be in m
	#serie must be ordered!
	sampling_spacing = x[1] - x[0]
	sampling_frequency = 1 / sampling_spacing
	nyquist_freq = 0.5*sampling_frequency
	print('Original nyquist frequency = ' + str(nyquist_freq) + 'cycles/m')
	
	n_samples_per_sub = np.fix(win_len * sampling_spacing) 
	
	start = x[0] + win_len / 2 * 100 #start/end positions in meters
	end = x[-1] - win_len / 2 * 100
	
	N_samples = np.fix((end - start) / samples_spacing)
	print('we are going to evaluate a total number of spectra: ' + str(N_samples))
	positions = np.linspace(start , end , N_samples) #positions of the sub-spectra
	print('Real spacing between spectra is ' + str(positions[1] - positions[0]))
	
	all_spectra = []
	for position in pb(positions):
		#get the subsignal:
		distances = np.abs(x - position)
		to_take = distances <= win_len / 2 #extract samples within the window (ids)
		subsignal = signal[to_take] #extract samples within the window (signal samples)
		p, f = mtspec.mtspec(subsignal,samples_spacing, pi, pad_to)
		p_to_keep = p[f <= maxfreq]
		p_to_keep = p_to_keep / np.max(p_to_keep)
		all_spectra.append(p)
	real_maxfreq = np.max(f[f<=maxfreq])
	all_spectra = np.array(all_spectra).T
	print('We have a final image with size: ' + str(all_spectra.shape))
	
		
	if plot == True:
		pylab.imshow(all_spectra, interpolation='nearest', aspect='auto', origin='Lower', extent=[start, end, 0, real_maxfreq])
		pylab.xlabel('position - cm')
		pylab.ylabel('frequency - cycles/m')
		pylab.colorbar()
		pylab.show()
	return ffts
	
	

def spectrogram_mtspec(pos, signal, win_len=1.2, pi=2, spacing=0.005, highfreq=8, pad_to=1000, sub_fact=10, plot=True):
	winlen_samples = np.int(win_len/spacing) #number of samples into a sample signal
	half_winlen = np.int(winlen_samples / 2)
	print('each window have size: ' + str(half_winlen * 2))
	ffts = []
	ns = []
	for n in pb(np.arange(0 + half_winlen, len(signal) - half_winlen)[::sub_fact]):
		ns.append(n)
		#working on signal
		small_signal = signal[n-half_winlen:n+half_winlen]
		#small_signal = mlab.detrend_linear(small_signal)
		#small_signal = small_signal * np.hanning(len(small_signal))
		#small_signal = np.concatenate((small_signal, np.zeros(pad_to-len(small_signal))))
		#getFFT
		import mtspec
		ps, f = mtspec.mtspec(small_signal, spacing, pi, pad_to)
		ps = ps / np.max(ps)
		
		ffts.append(ps[f<=highfreq])
	ffts = np.array(ffts)
	eval_pos = pos[ns]
	real_maxfreq = np.max(f[f<=highfreq])
	if plot == True:
		pylab.imshow(ffts.transpose(), interpolation='nearest', aspect='auto', origin='Lower', extent=[min(eval_pos), max(eval_pos), 0, real_maxfreq])
		pylab.xlabel('position - cm')
		pylab.ylabel('frequency - cycles/m')
		pylab.colorbar()
		pylab.show()
	return eval_pos, ffts


def mtm(x, K=3, NW=6, nFFT='default', plot=0, deltat=0.001, jkCIProb = 0.95, maxf=10, plotftest=True):
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
		ax1.plot(freqs, upperCI,'r:')
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


	
	