from __future__ import division
import numpy as np

from matplotlib.pyplot import plot, interactive, subplot, figure, title, xlabel, ylabel, grid, yscale
interactive(True)

from copy import deepcopy

# from ..NumericalMethods.SignalFilters import bandpass, spectrogram_mtspec, Spectrogram
#
# from scipy import signal
#
# import pylab


#import pdb; pdb.set_trace()

#from ..Elements import TimeSeriesXY

    

class TimeSeriesBasic(): 
    def __init__(self, y, title="series", unit = "m"):
        self.original_mean_ = y.mean()
        self.original_std_ = y.std()
        self.y_ = y
        self.title_ = title
        self.unit_ = unit

    def __repr__(self):
        return ( self.__module__ + ":\n" +  str(self.y_) )

        
    def remapValues(self, dictionary=dict()):
        """
        remaps the values of the y_  vector 
        the dict specify the mapping
        this is useful when dealing with discrete time-series. Like 
        lithologs in whcih the y_ vector could be somthing like
        {0,1,2,3} - four lithologies
        We may whant to put to zero also the lithology "1". 
        You can do that in this way:
        newseries = serie.remapValues({1:0})

        If dictioanry = {} (default) it returns deep-copied series
        
        Note: the series is deep-copied
        
        See Also: getUniqueValues()
        
        
        """
        copied = self.getDeepCopy()
        for k in dictionary:
            ids = np.argwhere(k==copied.y_)
            copied.y_[ids] = dictionary[k]
        return copied
            

    def printTitleInfos(self) :
        self.printInfosSeparator("Title")
        print("Series Name: \t" + str(self.title_) )
       
    def getNumberOfSamples(self):
        return self.y_.size
        
    def printInfosSeparator(self, title = "info", total_lenght = 100):
        print("\n")
        n_title = len(title)
        
        n = total_lenght - n_title - 2
        
        line_lenght = np.floor( n / 2 )
        line = '-' * line_lenght
        
        print(line  + " " +  title + " " + line)  

    def getDeepCopy(self):
        from copy import deepcopy
        return  deepcopy(self)
        
    def getMean(self):
        return self.y_.mean()
   
    def getStd(self):
        return self.y_.std()

    def __add__(self, other):
        assert(self.getNumberOfSamples() == other.getNumberOfSamples())
        from copy import deepcopy
        new = deepcopy(self)
        new.y_ = self.y_ + other.y_
        return new
        
    def getY(self):
        return self.y_

        
    def getUniqueValues(self):
        return np.unique(self.y_)
              
            
    def changeMean(self, mean=0.0):
        '''
        change the mean to the given value
        '''
        self.y_ = self.y_ - self.getMean() + mean

    def getMean(self):
        mdat = np.ma.masked_array(self.y_,np.isnan(self.y_))
        mm = np.mean(mdat)
        return mm

    def getStd(self):
        mdat = np.ma.masked_array(self.y_,np.isnan(self.y_))
        mm = np.std(mdat)
        return mm

    def getWithoutNans(self):
        from copy import deepcopy
        newseries = deepcopy(self)
        mean = newseries.getMean()
        newseries.y_ -= mean
        newseries.y_[np.isnan(newseries.y_)] = 0
        newseries.y_ += mean
        return newseries

    def changeStd(self, std= 1.0):
        '''
        change the std to the given value
        '''
        mean = self.getMean()
        self.y_ -= mean
        self.y_ = self.y_ / self.getStd() * std
        
        self.y_ += mean
        
    def restoreMeanStdToOriginal(self):
        self.changeMean(self.original_mean_)
        self.changeStd(self.original_std_)
        
    def normalizeMeanStd(self ):
        self.changeMean(0.0)
        self.changeStd(1.0)
        
    def getSeries(self):
        return self.y_
   
    def getDeepCopy(self):
        return deepcopy(self)

    def getPicker(self):
        figure()
        self.plot()
        from cyclopy.Orbitals import CollectPickPoints2
        p = CollectPickPoints2()
        return p

    def resampleAt(self, positions, step=None, method='linear', h=1):
        '''
        h is the h parameters for kernel smoother, unused for other methods
        methods:
        'linear'
        'rbf'
        'ks'

        positions is a array-type, with the positions at which resample the serie
        '''
        import scipy.interpolate

        x =self.getPositionVector()
        y = self.getY()

        if method == 'linear':
            interpolator = scipy.interpolate.interp1d(x, y, kind='linear', bounds_error=False, fill_value=np.nan)

        if method == 'rbf':
            interpolator = scipy.interpolate.Rbf(x, y, function='thin_plate')

        if method == 'ks':
            from cyclopy.NumericalMethods import KernelSmoother
            smoother = KernelSmoother.KernelSmoother(x, y)


        #call interpolators
        if method != 'ks':
            new_y = interpolator(positions)
        else:
            new_y = smoother(positions, h)

        from cyclopy.Elements import TimeSeriesXY, TimeSeriesEven



        new_series = TimeSeriesXY(positions, new_y, self.title_, self.unit_)

        return new_series


   
        