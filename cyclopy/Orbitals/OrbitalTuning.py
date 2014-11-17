# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 16:03:40 2014

@author: luca
"""

from matplotlib.pyplot import ginput, connect, gcf, axvline, disconnect, draw
import numpy as np
from scipy.interpolate import UnivariateSpline


class Picker:
    def __init__(self):
        self.finished = False
        self.picked = []
        self.lines = []
        self.connection = None

    def __call__(self, event):
        if event.inaxes and event.button is 1:
            clickX = event.xdata
            clickY = event.ydata
            current_fig = gcf()
            l = axvline(clickX)
            self.lines.append(l)
            self.picked.append((clickX, clickY))

            print("clicked on %f, %f" % (clickX, clickY))

        elif event.button is 2:
            self.close()
            print("middle button pressed. exiting")

        elif event.button is 3:
            if len(self.picked) == 0:
                print("nothing to undo")
                return

            self.picked.pop()
            self.lines.pop().remove()
            draw()

    def close(self):
        self.finished = True
        disconnect(self.connection)
        for l in self.lines:
            l.remove()

        draw()

    def asAgeModel(self, period = 100, zero = 0):
        x_pos = np.array(self.picked)[:,0]
        age_mod = AgeModelFromPicksAndPeriod(x_pos, period, zero=zero)
        return age_mod



def CollectPickPoints():
    print("Middle mouse button click to exit")
    data = np.array(ginput(0))[:, 0]
    return data


def CollectPickPoints2():
    p = Picker()
    connection = connect('button_press_event', p)
    p.connection = connection
    return p


def AgeModelFromPicksAndPeriod(picks, period=100, zero=0):
    ages = np.arange(len(picks)) * period + zero
    model = np.array([picks, ages]).T

    return model


def DepthToAge(depths, agemodel, k=1):
    interpolator = UnivariateSpline(agemodel[:, 0], agemodel[:, 1], k=k)  # k=1 means linear interpolation actually
    ages = interpolator(depths)
    return ages
    