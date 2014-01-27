# -*- coding: utf-8 -*-
"""
Created on Thu May 31 18:45:25 2012

@author: luca
"""

from __future__ import division
import SSA as SSA #import ssa methods
import numpy as np
from matplotlib.pyplot import axis, errorbar, close, plot, interactive, subplot, figure, title, xlabel, ylabel, interactive, grid, semilogy, hold, legend
interactive(True)

#load mratest.dat
mratest = np.loadtxt("mratest.dat")
x = mratest[:,0]

close('all') #close all

N = 300
#---------------------------------------------------------------
# Recreate Figure 1: the test signal alone and with noise added.
#---------------------------------------------------------------

figure()
plot(x, 'r', label="signal with noise")
hold (1)
plot(mratest[:,2], 'y', label="signal alone")
title('Figure 1: Test signal alone, and with added noise.')
legend()

#---------------------------------------------------------------
# Recreate Figure 3:
# This differs slightly from Allen and Smith's because I use the
# estimated AR(1) parameters rather than the known parameters.
#---------------------------------------------------------------

ssa = SSA.SSA(x)
ssa.setEmbeddingDimension(40)
ssa.setNumberOfMonteCarloSimulations(N)
ssa.update()
ssa.doMonteCarloSimulations()
c = ssa.getConfidence()

figure()
semilogy(ssa.Val_, '+r')
semilogy(c, 'c')
title('Figure 3: Eigenspectrum of test series and surrogate projections.')


#---------------------------------------------------------------
# Recreate Figure 4:
# This differs slightly from Allen and Smith's because I use the
# estimated AR(1) parameters rather than the known parameters.
#---------------------------------------------------------------


mE, fE = SSA.EOFDominantFrequencies(ssa.Eig_)

figure()
cm = np.mean(c, axis = 1)
semilogy (fE, ssa.Val_, '+r')

err = np.abs(c[:,0] - c[:,1]) / 2
errorbar (fE, cm, yerr=[c[:,1] - cm, cm - c[:,0]], linestyle='.')
axis([0, 0.5, 0.05, 15])
title("Fig 4")






