from cyclopy.NumericalMethods import KernelSmoother
from cyclopy.Elements.TimeSeriesBasic import TimeSeriesBasic

import numpy as np
from matplotlib.pyplot import title, xlabel, ylabel, plot, grid

from .TimeSeriesBasic import TimeSeriesBasic
from .TimeSeriesEven import TimeSeriesEven


class TimeSeriesXY(TimeSeriesBasic):
    """
    A time series made of two 1d vectors x, and y
    x are the sampling positions
    y the samples values
    args and kwargs are passed to TimeSeriesBasic. Please see docs fot that
    class.
    
    Note: The series is ALWAYS sorted so that x's values are increasing
    """

    def __init__(self, x, y, *args, **kwargs):
        assert (x.size == y.size)  # x and y must have same size
        #we want these two series to be sorted

        TimeSeriesBasic.__init__(self, y, *args, **kwargs)
        self.x_ = x

        self.sortSamples()  # ensure samples are sorted

    def plot(self):
        title(self.title_)
        xlabel(self.unit_ + " (" + str(self.getTotalRange()) + self.unit_ + " )")
        ylabel(self.title_)
        plot(self.x_, self.y_)
        grid()

    def addShift(self, shift):
        self.x_ += shift

    def getMinX(self):
        return min(self.x_)

    def getTotalRange(self):
        return max(self.x_) - min(self.x_)

    def sortSamples(self):
        sorted_ids = np.argsort(self.x_)
        self.x_ = self.x_[sorted_ids]
        self.y_ = self.y_[sorted_ids]

    def getMinMaxSamplingStep(self):
        s_steps = self.getSamplingSteps()
        return np.min(s_steps), np.max(s_steps)

    def getAverageSamplingStep(self):
        s_steps = self.getSamplingSteps()
        return s_steps.mean(), s_steps.std()

    def getSamplingSteps(self):
        n_steps = self.getNumberOfSamples() - 1
        s_steps = np.arange(n_steps, dtype=float)

        for i in np.arange(n_steps):
            s_steps[i] = self.x_[i + 1] - self.x_[i]
        return s_steps





    def printStepInfos(self):

        self.printInfosSeparator("Sampling Step")
        min_step, max_step = self.getMinMaxSamplingStep()
        avg_step, std_step = self.getAverageSamplingStep()

        print("Min Sampling Step: \t" + str(min_step))
        print("Max Sampling Step: \t" + str(max_step))
        print("Average Sampling Step: \t" + str(avg_step))
        print("Std of average Sampling Step:\t" + str(std_step))

    def printInfos(self):
        self.printTitleInfos()
        self.printStepInfos()

    def getAsEvenSpaced(self, step=None, method='linear', other=None):
        '''
        other is the h parameters for kernel smoother, unused for other methods
        
        NOTE this method should be rewritten to make use of the ResamplAt method, 
        that is a generic implementation
        '''
        import scipy.interpolate

        if step == None:
            step, std = self.getAverageSamplingStep()

        #compute the new x vector
        new_x = np.arange(np.min(self.x_), np.max(self.x_), step)

        if method == 'linear':
            interpolator = scipy.interpolate.interp1d(self.x_, self.y_, kind='linear')

        if method == 'rbf':
            interpolator = scipy.interpolate.Rbf(self.x_, self.y_, function='thin_plate')

        if method == 'ks':
            smoother = KernelSmoother(self.x_, self.y_)


        #call interpolators
        if method != 'ks':
            new_y = interpolator(new_x)
        else:
            new_y = smoother(new_x, other)

        new_series = TimeSeriesEven(new_y, step, np.min(self.x_), self.title_, self.unit_)

        return new_series

    def resampeAt(self, positions, step=None, method='linear', h=1):
        '''
        h is the h parameters for kernel smoother, unused for other methods
        methods:
        'linear'
        'rbf'
        'ks'

        positions is a array-type, with the positions at which resample the serie            
        '''
        import scipy.interpolate


        if method == 'linear':
            interpolator = scipy.interpolate.interp1d(self.x_, self.y_, kind='linear')

        if method == 'rbf':
            interpolator = scipy.interpolate.Rbf(self.x_, self.y_, function='thin_plate')

        if method == 'ks':
            smoother = KernelSmoother.KernelSmoother(self.x_, self.y_)


        #call interpolators
        if method != 'ks':
            new_y = interpolator(positions)
        else:
            new_y = smoother(positions, h)

        new_series = TimeSeriesXY(positions, new_y, self.title_, self.unit_)

        return new_series
        

    

    
            
        