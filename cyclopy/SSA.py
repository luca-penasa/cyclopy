# -*- coding: utf-8 -*-
"""
Created on Mon May 28 14:06:01 2012

@author: luca
"""
from __future__ import division
import numpy as np
from scipy.linalg import toeplitz
from scipy.signal import lfilter
import pylab

from matplotlib.pyplot import plot, interactive, subplot, figure, title, xlabel, ylabel, interactive, grid, semilogy
interactive(True)


import progressbar as pbar
pb = pbar.ProgressBar()

import MonteCarloMethods as mcm

__name__ = "SSA"
__doc__ = """
This code provide a partial reimplementation of the original matlab code from eric@gi.alaska.edu (Eric Breitenberger)
"""

class SSA:
    def __init__(self, y, emb_dim=100, acorr_method="unbiased"):
        """
        default init
        """
        self.emb_dim_ = None
        self.autocorrelation_method_ = "unbiased"
        self.num_montecarlo_ = 100
        
        #some checks on input signal
        assert(len(y.shape) == 1) # a vector?
        assert(len(y) != 0) # is void?
        
        self.y_ = y # time series
        self.n_samples_ = len(y) # number of samples
        
        self.setEmbeddingDimension(emb_dim)
        self.setAutocorrelationMethod(acorr_method)
        
        self.need_recomputation_ = True #always recompute everything when initializing
        
    def setEmbeddingDimension(self, emb_dim):
        """
        the embedding dimension to be used
        """
        assert((emb_dim -1) <= self.n_samples_) #not a too big dimension
        assert( type(emb_dim) == int ) #is it an int?
        
        if emb_dim != self.emb_dim_: #if emb_dim was changed we must recompute all
            self.emb_dim_ = emb_dim #set it
            self.need_recomputation_ = True #switch flag 
            
    def setAutocorrelationMethod(self, method):
        """
        methods to be used for estimating autocorrelation of your signal
        autocorrelation is then used to simulate a "similar" time-series of random values

        Accepted methods are:
        - "biased": Yule-Walker
        - "unbiased": traditional method
        see documentation for LagAutoCovariance
        """
        assert((method=="biased") or (method=="unbiased")) # only two methods implemented

        if self.autocorrelation_method_ != method: # is it the same?
            self.autocorrelation_method_ = method # set it
            self.need_recomputation_ = True # change flag
    
    def update(self):
        assert ((type(self.emb_dim_) == int) & (self.emb_dim_ > 0) ) # we need the embedding dimension to be set

        #do real computations
        self.Eig_, self.Val_ = EigenDecomposition(self.y_, self.emb_dim_, self.autocorrelation_method_) 
        self.Pcs_ = PrincipalComponents(self.y_, self.Eig_)
        self.Rcs_ = ReconstructedComponents(self.Pcs_, self.Eig_)
        self.need_recomputation_ = False #set the flag to False
        
    def getRcs(self, ids):
        """
        return the Reconstructed Components requested in ids
        you can use r_[1:20] for get all the RCS from 1 to 20
        or a list as [0,2,5,10] for a subset
        """
        if self.need_recomputation_ == True:
            self.update()
            
        return self.Rcs_[:,ids]
        
    def getSummedRcs(self, ids):
        """
        same as getRcs but return the complete reconstructed signal
        see getRcs for syntax hints
        """
        rcs = self.getRcs(ids)
        return np.sum(rcs, axis=1)
        
    def plotEigenValues(self):
        #is everything computed?
        if self.need_recomputation_ == True:
            self.update()
        
        figure()  
        semilogy( np.arange(self.emb_dim_) + 1, self.Val_, 'o')
        xlabel("Eigenvalue number")
        ylabel("Variance")
        grid()
        title("Eigenspectrum")
        
    def doMonteCarloSimulations(self):
        if self.need_recomputation_ == True:
            self.update()
        #parameters for the surrogate data
        gamma, alpha, mu2 = mcm.estimateAR1(self.y_)
        print ("Got estimates for surrogate data parameters:")
        print ("Gamma: " + str(gamma))
        print ("Alpha: " + str(alpha))
        
        L = np.zeros((self.emb_dim_, self.num_montecarlo_))
        vs = np.zeros(self.num_montecarlo_)
        s = np.zeros(self.num_montecarlo_)

        #produce surrogate data -> EACH ROW IS A REALIZATION
        self.surrogates_ = mcm.generateAR1Noise(self.n_samples_, self.num_montecarlo_, gamma, alpha)
        #self.surrogates_ = np.loadtxt('/home/luca/Desktop/Smirra/Code/randomX.txt')
	
        for i in pb(np.arange(self.num_montecarlo_ )):
            acv = LagAutoCovariance(self.surrogates_[:,i], self.emb_dim_ - 1, self.autocorrelation_method_)
            T = toeplitz(acv)
            #do projection of T on pre-computed eigenbasis -> USING FIXED BASIS!!
            V = self.Eig_.T.dot(T).dot(self.Eig_)
            
            s[i] = self.surrogates_[:,i].dot(self.surrogates_[:,i]) #variance of each surrogate
            vs[i] = np.sum(np.sum(V**2)) #squared Frobenius norm of V
            L[:, i] = np.diag(V) #get the diagonal and put in a column of L
        
        #put computed things where we can find them
        self.L_ = L
        self.s_ = s
        self.vs_ = vs
        
    def getConfidence(self, p = 0.95):
        #sorting L
        L_sorted = np.sort(self.L_, axis=1) #in each row we have the same lambda for all the realizations

        #now create the two columns vector, for upper and lover confidence
        #one for each eigenvalue
        c = np.zeros((self.emb_dim_, 2))
          
        #As in the original implementation by Breitenberger we do this:
        #% Now pull out the confidence limits: unless the limits
        #% fall exactly on an integer, this is done by
        #% selecting the values that lie closest to the *outside*
        #% of the selected confidence interval. For example, for
        #% N=100 and p=.95 (defaults) the 3rd and 98th sorted
        #% observation will be selected as the confidence limits.

        #probability at tails:
        p = (1-p) / 2 #two tailed
        #that corresponds to a "cut-off index", that depends on the number of realizations:
        p_index = self.num_montecarlo_ * p
        
        lower = np.floor(p_index) + 1
        upper = self.num_montecarlo_  - np.floor(p_index)
        
        c[:,0] = L_sorted[:, lower - 1] #-1 couse we have a zero index!
        c[:,1] = L_sorted[:, upper - 1]
          
              

          
        return c
        
    def setNumberOfMonteCarloSimulations(self, number_MC):  
        assert(number_MC > 30) #a simple limit
        self.num_montecarlo_ = number_MC
         

#TODO da sistemare questa non sono sicuro che sia corretta    
#    def computeInformationCriterions(self):
#        """
#        AIC and MDL
#        """
#        if self.need_recomputation_ == True:
#            self.update()
#            
#        p = len(self.Val_)
#        V = np.flipud(self.Val_)
#        L = np.zeros(p)
#        
#        nrm = np.flipud(np.arange(p))
#        sumlog = np.cumsum(np.log(V))
#        sumlog = np.flipud(sumlog)
#        logsum = np.log(np.cumsum(V))
#        logsum = np.flipud(logsum)
#        
#        L = self.n_samples_ * nrm * (( sumlog / nrm) - logsum + np.log(nrm))
#        
#        pen = np.arange(p)
#        pen = pen * (2*p - pen)
#        
#        aic = -L + pen
#        mdl = -L + pen * np.log(self.n_samples_) / 2
#        
#        kaic = np.arange(p)[aic == np.min(aic)]
#        kmdel = np.arange(p)[mdl == np.min(mdl)]
#        
#        
#        return aic, mdl
        
        
        
            
        

def LagAutoCovariance(y, n_lags, method = "unbiased"):   
    assert(n_lags < len(y))
    N = len(y) #number of samples in the series
    y = y - np.mean(y) #detrend
    #autocov vector will have lenght lag + 1 
    c = np.zeros(n_lags + 1) #create it
    
    
    for i in np.arange(n_lags + 1):
        c[i] = np.sum( y[0:N-i] * y[i:N] )
        
    if method == "biased": #Yule-Walker
        c = c / N
        
    elif method == "unbiased":
        c = c / ( - np.arange(n_lags + 1) + N)
        
    return c
        
    
def EigenDecomposition(y, M, method="unbiased"):
    """
    method is passed to LagAutoCovariance for the lag-covariance computation
    """
    #get autocov
    c = LagAutoCovariance(y, M - 1, method)
    #as Toeplitz matrix:
    T = toeplitz(c)
    #now perform Eigen decomposition
    V, E = np.linalg.eig(T)
    # V -> eigenvalues
    # E -> eigenvectors
    
    #now sort eigenvalues in decreasing order
    sort_id = np.argsort(V) #ids sorting eigenvalues in creasing order
    sort_id = np.flipud(sort_id) #same in decreasing
    Val = V[sort_id] #perform real sorting
    Eig = E[:, sort_id] #also for E
    
    return Eig, Val
    
    
def PrincipalComponents(y, Eig):
    '''
    each row is a PC
    E: eigenvectors
    '''
    N = len(y) #number of samples
    #zero-mean
    y = y - np.mean(y)
    M, tmp = Eig.shape
    assert (M == tmp) #must be square
    
    PCs = np.zeros((N-M+1, M))
    #now reconstruct each PC
    for i in np.arange(N-M+1): 
        w = y[i:i+M]
        PCs[i, :] = w.dot(Eig)
    
    return PCs
    
def ReconstructedComponents(PCs, Eig):
    M, tmp = Eig.shape    
    pc_row, pc_col = PCs.shape
    
    #checks:
    assert(M == tmp)    
    assert(pc_col == M)
    
    N = pc_row + M - 1 # Assumes A has N-M+1 rows.
    
    R = np.zeros((N, M))
    Z = np.zeros((M-1, M))
    
    PCs = np.vstack([PCs, Z])
    
    #now compute the RCs
    for k in np.arange(M):
        R[:,k] = lfilter(Eig[:,k], M, PCs[:,k])
       
    #Adjust first M-1 rows and last M-1 rows
    for i in np.arange(M-1):
        R[i, :] = R[i, :] * (M/(i+1)) 
        R[N-i-1, :] = R[N-i-1, :] * (M/(i+1))
        
    return R
        
 
def EOFDominantFrequencies(Eig, n_max = 500):
    """
    Find the dominant frequency for each EOF
    Eig -> Matrix of Eigenvectors 
    n_max -> max number of frequencies to be checked   
    
    Uses a Reduced Fourier Transform for computing coefficients at requested freqs
    """
    s = EOFSymmetry(Eig)
    M, K = Eig.shape
    
    #de-mean eigenfunctions
    Eig = Eig - np.ones((M , 1)) * np.mean(Eig, axis=0) 
    
    f = np.array([np.linspace(0, 0.5, n_max)]).T #freqencies to be checked
    F = np.zeros((n_max, K)) #values at that freq - to be filled in
    
    j = (np.arange(M) + 1) - (M + 1) /2 
    Cc = f * j
    Cs = np.sin(2 * np.pi * Cc)
    Cc = np.cos(2 * np.pi * Cc)
    
    r2c = np.sum( (Cc ** 2).T , axis=0) / M
    r2s = 1 - r2c
    
    r2s[0] = np.spacing(1) #no zero divide in this way!
    
    if r2s[n_max - 1] == 0:
        r2s[n_max - 1] == np.spacing(1)
    
    
    for k in np.arange(K):
        if s[k] == 1:
            F[:,k] = np.abs(Cc.dot(Eig[:,k])**2 / r2c)
            
        elif s[k] == 0:
            F[:,k] = np.abs(Cs.dot(Eig[:,k])**2 / r2s)
        
        elif s[k] == -1:
            F[:,k] = np.abs(Cc.dot(Eig[:,k]))**2 + np.abs(Cs.dot(Eig[:,k]))**2
            
    
    mrft = np.max(F, axis=0) / M
    f = np.argmax(F, axis=0) #da cambiare -> devo trovare il max in ogni colonna
    
    
    f = 1/2 * (f) / (n_max-1)
    
    return mrft, f
    
        
    
    
    
    
    
    

def EOFSymmetry(Eig, tol = 10**4*np.spacing(1)):
    """
    Check the simmetry of the Empirical Orthogonal Functions
    A vector of this type is returned, as in the Breitenberger implementation:
    % If the i-th EOF is symmetric,  s(i)=1,
    % if anti-symmetric,             s(i)=0.
    % if neither sym. or anti-sym., s(i)=-1.
    %
    % s(i)=-1 is only possible if a non-Toeplitz (BK type)
    % covariance matrix was used, or if the tolerance 'tol'
    % is not set high enough. 'tol' is set to tol=10^4*eps
    % by default, or it can be specified as the second argument.
    %
    % If 'tol' is set to the string '1 or 0' or ('0 or 1') the
    % output will be forced to give only ones and zeros. EOFSYM
    % will decide whether the EOFs are symmetric or anti-symmetric
    % based on which assumption gives the lowest rms error.
    """
    M, K = Eig.shape
    s = np.zeros(K) #this will be the output
    
    if type(tol) != str: #three cases are possible
        for k in np.arange(K):
            if M%2 == 0:
                L = Eig[0:M/2, k]
                R = np.flipud(Eig[M/2: M, k])
            else:
                L = Eig[0:(M-1)/2 , k]
                R = np.flipud(Eig[M/2 + 1: M, k])
                
            if np.max(np.abs(L-R)) < tol:
                s[k] = 1
            elif np.max(np.abs(L+R)) < tol:
                s[k] = 0
            else:
                s[k] = -1
                print("EOF " + str(k) + " is neither symmetric or antisymmetric")
                
    elif (tol == "1 or 0") or (tol == "0 or 1"):
        Esym = np.sum(np.sqrt( np.flipud(Eig) - E  ) ** 2 )
        Easym = np.sum(np.sqrt( - np.flipud(Eig) - E  ) ** 2 )
        s = (np.sign/(Easym - Esym) + 1) / 2
        ties = pylab.find(s == 1/2) #resolve any ties -> call them symmetric!
        s[ties] = np.ones(len(ties))
    else:
        print("Improper specification of the tolerance")
        
        
    return s            
    

    
    
