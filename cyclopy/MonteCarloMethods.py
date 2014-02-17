# -*- coding: utf-8 -*-
"""
Created on Mon May 28 20:01:05 2012

@author: luca
"""
from __future__ import division
import numpy as np
from scipy.optimize import fsolve as fsolve
from scipy.optimize import bisect as bisect
from numpy.random import randn



__name__="MonteCarloMethods"


def estimateAR1(x):
    """
    Given an AR(1) realization estimates the parameters gamma, alpha and mu2
    for the model.
    """
    c0 = estimateAR1LagCov(x, 0)
    c1 = estimateAR1LagCov(x, 1)
    
    g0 = c1/c0 #intial estimate of gamma
    N = len(x)    
    
    #now usa a nonlinear solver for finding the zero of gammaEstimatorFunction
    
    gamma_est = fsolve(gammaEstimatorFunction, g0, args=(N, g0), xtol=1e-10)
    
    #compute the mu2
    N_part = N - (np.arange(N-1) + 1)
    gamma_part = gamma_est**(np.arange(N-1) + 1)
    
    mu2_est = 1/N + 2/(N**2) * np.sum(N_part*gamma_part)
    
    #and get the estimated c0 and alpha for return:
    c0_est = c0 / (1-mu2_est)
    alpha_est = np.sqrt((1 - gamma_est**2) * c0_est)
    
    return gamma_est[0], alpha_est[0], mu2_est
    

def gammaEstimatorFunction(gamma, N, g0):
    """
    N lenght of the series
    gamma the a gamma value
    g0 is the ratio c1/c0, that is also the inital estimate of gamma
    This function implements equations 9 and 11 from Allen & Smith 1996
    These are used for the gamma estimation of the AR(1) noise model
    Root of this function must be find with some nonlinear method
    Return is the current value of the function
    """
    
    #compute mu_squared for this gamma value
    N_part = N - (np.arange(N-1) + 1)
    gamma_part = gamma**(np.arange(N-1) + 1)
    
    mu2 = 1/N + 2/(N**2) * np.sum(N_part*gamma_part)
    
    function_value = g0 + mu2 * (1 - g0) - gamma
    
    print ("Iteration for finding gamma: " + str(gamma) + " with function value: " + str(function_value))
    
    return function_value
    
    

def estimateAR1LagCov(x, lag):
    """
    Estimate the Lag-l covariance with a natural estimator
    for any given time series x
    See eq 8 of Allen & Smith 1996
    This estimator is less biased than Yule-Walker estimate,
    but it is still biased 
    """
    x = x - np.mean(x) #remove mean
    N = len(x)
    
    assert(lag >= 0)
    
    
    if lag == 0:
        a = 1/N
        b = np.sum(x*x)

    if lag > 0:
        a = 1/(N-lag)
        b = np.sum( x[0:N-lag] * x[lag:N] )
        
    
    c_lag = a*b
    return c_lag
    

def generateAR1Noise(n, c, gamma, alpha):
    """
    generate red noise
    from a AR(1) system, 
    n number of samples
    c number of time series to be generated
    gamma lag-1 autocorrelation
    alpha noise innovation variance parameter
    
    return an n x c matrix of red noise series
    """
    X = np.zeros((n, c))
    X[0, :] = np.sqrt(alpha**2 / (1 - gamma**2)) * randn(c)
    z = alpha * randn(n, c)
    
    for i in np.r_[1:n]:
        X[i, :] = gamma * X[i - 1, :] + z[i, :]
        
    X = X - np.ones((n, 1)) * np.mean(X)
    
    return X

    
    
    
    
