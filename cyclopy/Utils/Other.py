# -*- coding: utf-8 -*-
"""
Created on Tue May 15 11:01:14 2012

@author: luca
"""
import numpy as np
from matplotlib.pyplot import figure, plot

from copy import deepcopy

import pylab

from cyclopy.Elements import TimeSeriesXY

class TiePoints(): #Base class  for each tie point type
    def __init__(self, x, y):
        '''
        x -> stratigraphic positions or time of the original time series
        y -> a variable for that tie point
        
        in this latter case we use the frequency domain information for estimating a sedimentation rate and how it varies along the series
        in the first case we use the x to y constrains for estimating a model for streching/compressing the series
        '''
        self.x_ = x
        self.y_ = y
        
    def __repr__ (self):        
        string = "x: " + str(self.x_) + "\n" + "y: " +  str(self.y_)
        return string
        
      
     
    def plot(self, newfig=True):
        if newfig:
            figure()
        
        plot(self.x_, self.y_, '-o')
            
    def getAsTimeSeriesXY(self):
        return TimeSeriesXY(self.x_, self.y_)
        
    def getTiesInBetween(self, limit ):
        '''
        permits to get a copy of the ties, but only the ones between two limits
        limit is a 2 element vector or list [min_x, max_x]
        '''
        copied = deepcopy(self)
        ids = pylab.find((self.x_ > limit[0]) & (self.x_ < limit[1]))
        assert (len(ids) != 0)
        copied.x_ = self.x_[ids]
        copied.y_ = self.y_[ids]
        return copied

      

class FrequencyTiePoints(TiePoints):
    def __init__ (self, x, y):
        '''
        Each tie point is interpreted as a position (x) and a frequency (y) of the original time series
        '''
        TiePoints.__init__(self, x, y)
        

         
    def setRequestedFrequency(self, freq=1):
        '''
        the frequency at which we need to tune the observed frequencies contained into y
        '''
        self.freq_ = freq
    
    def getSedRateModel(self):
        rate = self.freq_ / self.y_ #sedimentation rate for each tie
        

       

class SimpleModel():
    def __init__(self, poly = np.array([0,0])):
        '''
        an age model gives a value for any position (and vice-versa)
        it is modeled as a polynomial with polyval      
        Using more than one model we can create a composite model made up of several intervals,
        each one indipendently estimated in some way
        '''
        #default values
        self.p_ = poly
        self.d_ = 1
        self.bounds_ = [0,0]
        
        have_changed_ = True
        
    def getValue(self, position):
        '''
        ask to the model the age of a given stratigraphic position
        '''
        if (self.have_changed == True):
            self.update()
            
        return np.polyval(self.p_, position)
    
    def getPosition(self, age):
        '''
        ask to the model where is a given age
        this can be done only for 1-degree models
        that are not "flat"
        '''
        if (self.have_changed == True):
            self.update()
            
        assert (self.d_ == 1)
        assert (self.p_[0] != 0.0)
        return (age - self.p_[1]) / self.p_[0]
        
    def setModelDegree(self, degree):
        '''
        for changing the degree of the poly
        '''
        self.d_ = degree
        self.have_changed = True
    
    
    def setTies(self, ties):
        self.ties_ = ties
        #update bounds for these ties
        self.bounds_ = [np.min(self.ties_.x_), np.max(self.ties_.x_)]
        self.have_changed = True
        
    def update(self):
        '''
        for fitting this linear model we need a series of known Tie points 
        points for which the modeled variable is know
        returns squared residuals of fitting
        '''
        if (self.have_changed == True):
            self.p_ = np.polyfit(self.ties_.x_, self.ties_.y_, self.d_)
            
        self.have_changed = False
            
    
    def getResiduals(self):
        
        return self.getValue(self.ties_.x_) - self.ties_.y_
        
    def plot(self, newplot=True):
        if self.have_changed == True:
            self.update()
            
        if newplot == True:
            figure()
            
        a = np.array(self.bounds_)
        b = self.getValue(a)
        plot(a, b)
        
        
class CompositeModel():
    def __init__(self, ties = None, breaks=None):
        '''
        init permits to set ties and breaks all at once or with the dedicated methods
        ties is a TiePoints-derived Object
        breaks are just a list or numpy array with the positions of breaks
        '''
        self.ties_ = ties
        self.breaks_ = breaks        
        self.have_changed_ = True
        self.d_ = 1 #degree of poly
        
                
    def setTies(self, ties):
        self.ties_ = ties
        self.have_changed_ = True
    
    def setBreakPoints(self, breaks):
        self.breaks_ = breaks 
        #and sort them
        #from lower to higher
        self.breaks_ = np.array(breaks)[np.argsort (breaks)].tolist()
        self.have_changed_ = True
        
    def update(self):
        '''
        we separate the points depending the breakpoints
        and for each segment a SimpleModel is fitted
        '''
        #check some things
        assert (self.ties_ != None)
        assert (self.breaks_ != None)
        
        if (self.have_changed_ == False):
            return
            
        self.models_ = [] #delete the models
        
        #fit the first model
        #get pertinent ties
        tmp_ties = self.ties_.getTiesInBetween([-np.inf, self.breaks_[0]])
        
        model = SimpleModel()
        model.setTies(tmp_ties)
        model.setModelDegree(self.d_)
        model.update() #update it
        self.models_.append(model)
        
        for m in np.r_[1:len(self.breaks_)]: 
            model = SimpleModel()
            tmp_ties = self.ties_.getTiesInBetween([self.breaks_[m-1], self.breaks_[m]])
            model.setTies(tmp_ties)
            model.setModelDegree(self.d_)
            model.update()
            self.models_.append(model)
            
        #now the last model
        tmp_ties = self.ties_.getTiesInBetween([self.breaks_[len(self.breaks_) - 1], np.inf])
        model = SimpleModel()
        model.setTies(tmp_ties)
        model.setModelDegree(self.d_)
        model.update() #update it
        self.models_.append(model)
    
        self.have_changed_ = True

    def _split_(self, series, bounds):
        '''
        split a time series extracting only the points within the bounds
        return the indices of elements within the limits
        series -> a 1d array
        bounds -> list or array [min, max]
        '''
        return pylab.find((series > bounds[0]) & (series <= bounds[1]))
        
    
    def getValue(self, positions):
        N = len(positions) #number of points
        out = np.zeros(N)
        #for each point we want to evaluate we need to know the pertinent breakpoint 
        #first piece
        ids = self._split_(positions,[-np.inf, self.breaks_[0]])
        out[ids] = self.models_[0].getValue(positions[ids])
        
        #now do a for for the pieces in between
        for m in np.r_[1:len(self.breaks_)]: 
             ids = self._split_(positions,[self.breaks_[m-1], self.breaks_[m]])
             out[ids] = self.models_[m].getValue(positions[ids])
             
        #for the last piece
        ids = self._split_(positions,[self.breaks_[len(self.breaks_) - 1], np.inf])
        out[ids] = self.models_[len(self.breaks_) ].getValue(positions[ids])
        
        return out
        
    def getModelBounds(self):
        tmp = []
        for m in self.models_:
            tmp.append(m.bounds_)
            
        tmp = np.array(tmp)
        return [np.min(tmp), np.max(tmp)] 
    
    def plot(self, newfig = True, ties=True, model=True):
        '''
        do a plot of the model, the ties or both, within or not the current figure
        '''
        if newfig:
            figure()
        
        if ties:
            self.ties_.plot(False)
            
        if model:
            #then plot our model
            min, max = self.getModelBounds()
            tmp_x = np.linspace(min, max, 1000)
            this = self.getValue(tmp_x)
            
            plot(tmp_x, this)

            
        
