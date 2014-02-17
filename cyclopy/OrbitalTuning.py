# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 16:03:40 2014

@author: luca
"""

try:
    from matplotlib.pyplot import ginput
except ImportError:
    print("Need matplotlib to work. Install it - using pip install matplotlib for example")    
    
try:    
    import numpy as np
except ImportError:
    print("Need numpy to work!")    
    
    
try:
    from scipy.interpolate import UnivariateSpline
except ImportError:
    print("You also need scipy")            
   

#def extrap1d(interpolator):
#    xs = interpolator.x
#    ys = interpolator.y
#
#    def pointwise(x):
#        if x < xs[0]:
#            return ys[0]+(x-xs[0])*(ys[1]-ys[0])/(xs[1]-xs[0])
#        elif x > xs[-1]:
#            return ys[-1]+(x-xs[-1])*(ys[-1]-ys[-2])/(xs[-1]-xs[-2])
#        else:
#            return interpolator(x)
#
#    def ufunclike(xs):
#        return np.array(map(pointwise, np.array(xs)))
#
#    return ufunclike     
    
def CollectPickPoints():
    print("Middle mouse button click to exit")
    data = np.array(ginput(0))[:,0]
    return data
    
def AgeModelFromPicksAndPeriod(picks, period=100, zero=0):
    ages = np.arange(len(picks)) * period + zero
    model = np.array([picks, ages]).T
    
    return model
    
def DepthToAge(depths, agemodel):
    interpolator = UnivariateSpline(agemodel[:,0], agemodel[:,1], k=1) # k=1 means linear interpolation actually
    ages =  interpolator(depths)                          
    return ages
    