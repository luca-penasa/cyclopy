__author__ = 'luca'

import os



def getDataPath():
    datafolder = os.path.dirname(os.path.realpath(__file__))
    return datafolder

def getLa2010DataPath():
    folder = getDataPath()

    path = folder + os.path.sep +  'La2010'
    return path