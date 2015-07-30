__author__ = 'luca'

import numpy as np
import os

from cyclopy.Elements import TimeSeriesEven

from cyclopy.Data import getLa2010DataPath

datafolder = "../Data/La2010/"

def loadLa2010Eccentricity():
    file = getLa2010DataPath() + os.path.sep + "La2010a_ecc3L.dat"
    x,y = np.loadtxt(file).T

    ts = TimeSeriesEven(y[::-1], x_start=np.min(x), x_step=np.abs(x[1] - x[0]))
    ts.unit_ = "kyr"


    return ts
