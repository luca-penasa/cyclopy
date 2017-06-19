from cyclopy.Elements import TimeSeriesEven
import numpy as np

class TimeSeriesMerger():
    def __init__(self):
        self.series = [] #list of series
        self.output = None

    def addTimeSeriesToMerge(self, ts):
        self.series.append(ts)

    def initOutput(self):
        minx = +np.inf
        maxx = -np.inf

        for ts in self.series:
            max=ts.GetMinX()
            mix = ts.GetMaxX()