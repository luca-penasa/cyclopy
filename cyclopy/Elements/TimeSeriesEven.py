from cmath import sqrt
import scipy
import scipy.stats

from matplotlib.pyplot import title, figure, xlabel, ylabel, plot, grid, yscale

import numpy as np

from copy import deepcopy

from ..NumericalMethods import KernelSmoother

from ..NumericalMethods.SignalFilters import bandpass

from ..NumericalMethods.Spectrogram import Spectrogram

from ..NumericalMethods import SignalFilters


from ..Elements import TimeSeriesBasic


from scipy import signal

class TimeSeriesEven(TimeSeriesBasic):
    """
    A time series made of only a vector.
    optionally the samples start position and the sampling step may be 
    provided.
    aka an evenly spaced time series
    
    args and kwargs may be the same of a TimeSeriesBasic
    
    Note: Given most spectral methods are for evenly spaced ts this class
    support nice and advanced analysis
    
    See Also: TimeSeriesXY and TimeSeriesBasic
    """
    def __init__(self, y, x_step = 1.0, x_start=0.0, *args, **kwargs):
        TimeSeriesBasic.__init__(self, y, *args, **kwargs)
        self.x_step_ = x_step
        self.x_start_ = x_start

        
    def plot(self, nfig=False, *args, **kwargs):
        x = self.getPositionVector()        

        if nfig:
            figure()
        
        title(self.title_)
        xlabel(self.unit_)
        ylabel(self.title_)
        plot(x, self.y_, *args, **kwargs)
        grid()
        
    def addShift(self, shift):
        self.x_start_ += shift
        
    def getMinX(self):
        return self.x_start_

    def sampleAt(self, position, max_distance):
        x = self.getPositionVector()
        sq_distances = (x - position)**2
        argsorting = np.argsort(sq_distances)
        best_id =  argsorting[0]

        if sq_distances[best_id] <= max_distance**2:
            return self.y_[best_id]
        else:
            return np.nan


        
    def getMaxX(self):
        return np.max(self.getPositionVector())
        
    def TuneIt(self, agemodel, spline_k = 1):
        try:
            from ..Orbitals import DepthToAge
        except ImportError:
            print ("cannot use this method. cannot locate orbital tunig module")
        
        ages = DepthToAge(self.getPositionVector(), agemodel, k=spline_k)

        from ..Elements import TimeSeriesXY
            
        return TimeSeriesXY(ages, self.getY())

    def toTime(self, sedrate):
        """
        from space to time-domain
        sedrate in m/Ma
        new x-scale will be in kyr
        """

        srk=sedrate/1000; # converts sedrate to m/kyr
        self.x_step_/= srk;
        self.x_start_ /= srk
        self.unit_ = "kyr"

    def toSpace(self, sedrate):
        srk=sedrate/1000; # converts sedrate to m/kyr
        self.x_step_*= srk;
        self.x_start_ *= srk
        self.unit_ = "m"
        
    def crossCorrelate(self, series_b, do_plot=True, normalize=True):
        assert(self.x_step_ == series_b.x_step_)
        
        
        ccorr = np.correlate( self.y_, series_b.y_, mode='full')

#       N = max(self.y_.size, series_b.y_.size) - min(self.y_.size, series_b.y_.size) + 1        
#       N = max(self.y_.size, series_b.y_.size) 
        N = self.y_.size + series_b.y_.size -1
        initial_shift = self.x_start_ - series_b.x_start_

        shifts = (np.arange(N) - len(series_b.y_) ) * self.x_step_ + initial_shift        
        
        if do_plot:
            title("Cross correltation of series " + self.title_ + " vs. " + series_b.title_)
            plot(shifts,  ccorr)
            
            grid()
            
        best_shift = shifts[ccorr == np.max(ccorr)]
        self.printInfosSeparator("Cross correlation report")
        print("Shift giving best correlation: " + str(best_shift))
        
        return shifts, ccorr, best_shift
        
    def crossPlot(self, series_b, shift = 0.0):
        self.plot()
        shifted = deepcopy(series_b)
        shifted.x_start_ = shifted.x_start_ + shift
        shifted.plot()
        
    def bandpassFFT(self, fa, fb):
        new_series = deepcopy(self)
        new_series.y_ = bandpass(new_series.y_, fa, fb, new_series.x_step_)
        return new_series
        
    def bandpassZeroShift(self, fa, fb, n=201, show_res = False):
        f = SignalFilters.designLinearBandpass(fa, fb, self.x_step_, n, show= show_res)
        new_series = deepcopy(self)
        filtered = signal.filtfilt(f, [1.0], self.y_)
        
        new_series.y_ = filtered
        
        return new_series
        

    def getMTMSpectrum(self, pi = 2, pad_to=100000, norm=True, 
                       from_x=None, to_x=None, 
                       max_f=None, 
                       detrend=True, det_type='linear', 
                       stats = False,
                       p_crit = 0.95
                       ):
        """
        Get the MTM (Thompson multitaper method) spectrum for the series
        
        This is an interface to mtspec python module
        
        Keyword arguments:
        pi          --  the bandwhidth passed to mtspec.mtspec method        
        pad_to      --  if needed pad the series with zeros
        norm        --  normalize forcing the higher peak to 1.0
        from/to_x   --  compute the spectrum only for a sub-slice of the series
        max_f       --  get rid of frequencies above this value (cycles/time)
        stats       --  gives back a lot of additional stats and data from the 
                        mtspec routines
        p_crit      --  critical p-value for the f-test used for reshaped spectra
        """     
        
            
        import mtspec
        
        if ((from_x!=None) and (to_x!=None)):
            assert(from_x < to_x)
            ser = self.getSlice(from_x, to_x)
        else:
            ser = self
        
        
        if len(ser.y_) < pad_to:
            zeros = np.zeros(len(ser.y_ - pad_to))
            new_signal = np.concatenate([ser.y_ , zeros])
        else:
            new_signal = ser.y_
            
        if detrend:
            new_signal = signal.detrend(new_signal, type=det_type)
            

        ampl, freqs = mtspec.mtspec(new_signal, ser.x_step_, pi, nfft=pad_to)
        
        statistics = {}
        
        if stats:
            ampl_reshaped, freqs2, eigensp, eigencoeff, weights, jack, f, dof = mtspec.mtspec(new_signal, ser.x_step_, pi, statistics=True, fcrit=p_crit, optional_output=True, rshape=0, nfft=pad_to)
            nd = eigensp.shape[1]
            lev = scipy.stats.f.ppf(np.array([0.80,0.90,0.95,0.99]), 2, 2*(nd - 1))
            statistics = {      'reshaped': ampl_reshaped, 
                                'eigenspectra': eigensp,
                                'eigencoeff': eigencoeff,
                                'weights': weights,
                                'jack5_95': jack,
                                'f_values': f,
                                'dof': dof}
                                                        
        if max_f != None:
            ids = np.argwhere(freqs <= max_f)            
            ampl = ampl[ids]
            freqs = freqs[ids]
            
            if stats: # cut away also statistics
                for k in statistics.keys():
                    arr = statistics[k]
                    statistics[k] = arr[ids]                    
                
        
        if stats:
            statistics['p_levels'] = lev
                
        if norm == True:
            ampl = ampl/np.max(ampl)
                     
        
        if stats:
            return freqs, ampl, statistics
        else:
            return freqs, ampl
        
    def plotMTMSpectrum(self, nfig=False, ylog=False, **kwargs):
        """
        kwargs are passed to getMTMSpectrum. See that method for ref
        
        Keyword arguments
        nfig        --  Create  new figure before plotting (def=False)
        ylog        --  Use log10 scale for Y axis (def=False)
        
        Note: additional kwargs will be passed to getMTMSpectrum()
        """
        freqs, power = self.getMTMSpectrum(**kwargs)
        
        if nfig:
            figure()

        plot(freqs, power)
                
        if ylog:
            yscale('log')

        return freqs, power
                
        
    

            
    def getSetOfEvolutionarySpectra(self, positions, winsize=5.0, **kwargs):  
        spectra = []
        for  pos in positions:
            this_slice = self.getSlice(pos - winsize*0.5, pos+winsize*0.5)
            # figure()
            # this_slice.plot()
            f, sp = this_slice.getMTMSpectrum(**kwargs)
            spectra.append([f, sp])

            
            
        return spectra
        
    def plotEvolutionarySpectraLines3D(self, positions, mode='lines', **kwargs):
        """
        mode        -- possible are 'lines' or 'surface'
        """
        spectra = self.getSetOfEvolutionarySpectra(positions, **kwargs)
        import matplotlib.pyplot as plt

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        if mode=='lines':        
            i = 0
            for sp in spectra:    
                fre, power = sp
                pos = positions[i] * np.ones(len(fre)) 
                ax.plot(pos, fre, '-k', zs=power)            
                i += 1
        elif mode=='surface':
            x = np.array([])
            y = np.array([])
            z = np.array([])
            
            i = 0                            
            for sp in spectra: 
                fre, power = sp
                pos = positions[i] * np.ones(len(fre)) 
                x= np.concatenate([x, pos])                            
                y= np.concatenate([y, fre])
                z= np.concatenate([z, power])
                i += 1
            ax.plot_surface(x, y, z)
            
        
    def plotHilbertTransformAbs(self):        
        hil = self.getHilbertTransformAbs()
        plot(np.arange(len(self.y_)) * self.x_step_ + self.x_start_, abs(hil.y_) )
        
    def getHilbertTransformAbs(self):
        hil = deepcopy(self)        
        import scipy.signal
        hil.y_ = abs(scipy.signal.hilbert(self.y_))
        return hil
        
    def getNumberOfSamples(self):
        return len(self.y_)
    def getOverlappingInterval(self, series_b):
        Na = len(self.y_)
        Nb = len(series_b.y_)
        
        maxa = Na * self.x_step_ + self.x_start_
        maxb = Nb * series_b.x_step_ + series_b.x_start_
        
        mina = self.x_start_
        minb = series_b.x_start_
        
        upper = np.min([maxa, maxb])
        lower = np.max([mina, minb])
        
        if upper < lower:
            print('Series do not intersect each other')
            return 0
        else:
            return [lower, upper]
    
    def getPositionVector(self):
        return np.arange(len(self.y_)) * self.x_step_ + self.x_start_
        
    
    def getSlice(self, lower, upper):
        new_series = deepcopy(self)
        x = self.getPositionVector()
        all_ids = np.arange(len(self.y_))
        
        ids = all_ids[(x >= lower ) & (x <= upper) ]
        new_series.y_ = self.y_[ids]
        new_series.x_start_ = x[ids][0]
        
        return new_series
        
    def doCorrelationTest(self, series_b, method='pearson'):
        a, b = self.getOverlappingInterval(series_b)
        a_part = self.getSlice(a,b)
        b_part = series_b.getSlice(a,b)
        
        #the two MUST have the same size
        a_n = a_part.getNumberOfSamples()
        b_n = b_part.getNumberOfSamples()

	    #TODO finish this, using scipy.stats.pearsonr

    def getKS(self, h=1.0):
        smoother = KernelSmoother(self.getPositionVector() , self.y_)
        new_y = smoother(self.getPositionVector(), h)
        new_series = deepcopy(self)
        new_series.y_ = new_y
        return new_series
       
    def detrendKS(self, h=1.0):
        smoother = KernelSmoother(self.getPositionVector() , self.y_)
        new_y = smoother(self.getPositionVector(), h)
        new_series = deepcopy(self)
        new_series.y_ = new_series.y_ - new_y
        return new_series
        
    def plotCWT(self, order=6):
        from mlpy.wavelet import cwt, autoscales, fourier_from_scales
        import matplotlib.pyplot as plt
        omega0 = order
        #compute scales 
        scales = autoscales(N=self.getNumberOfSamples(), dt=self.getSamplingStep(), dj=0.05, wf='morlet', p=omega0)        
        
        spec = cwt(self.getSeries(), dt=self.getSamplingStep(), scales=scales, wf='morlet', p=omega0)
        
        freq = fourier_from_scales(scales, 'morlet', omega0)
        
        # approximate scales through frequencies
        #freq = (omega0 + np.sqrt(2.0 + omega0**2))/(4*np.pi * scale[1:])
        #get as fourier freqs
        
        
        plt.figure()
        ax1 = plt.subplot(211)
        ax2 = plt.subplot(212, sharex=ax1)
        
        t = np.arange( self.getNumberOfSamples() ) * self.x_step_ + self.x_start_
        
        ax1.plot(self.getPositionVector(), self.y_, 'k')
        extent=[t[0], t[-1], freq[0], freq[-1]]
        img = ax2.imshow(np.abs(spec), extent=extent, aspect='auto') 
        
        #plt.colorbar(img)
        plt.show()
        return extent,spec
        
    def getSpectrogram(self, win_len=1, pi=1, highfreq=1, detrend='linear', method='mtspec', pad_to=1000, sub_fact=1, plot=True, max_normalize=True):
        s = Spectrogram(self)
        s.setHighFreq(highfreq)
        s.setLocalDetrend(detrend)
        s.setMethod(method)
        s.setSubsamplingFactor(sub_fact)
        s.setUseMaxNormalization(max_normalize)
        s.setWindowLenght(win_len)
        s.setZeroPadding(pad_to)
        s.setMTMpi(pi)
        s.update()
        if plot:
            s.plot(True)
        return s
                
        
        eval_pos, ffts = spectrogram_mtspec(self.getPositionVector(), self.y_, win_len=win_len, pi=pi, spacing = self.x_step_, highfreq = highfreq, pad_to = pad_to, sub_fact = sub_fact, plot = plot, max_normalize = max_normalize, use_log = use_log)
        return eval_pos,ffts
    
    def getSamplingStep(self):
        return self.x_step_
    
    def setSamplingStep(self, step=1):
        self.x_step_ = step