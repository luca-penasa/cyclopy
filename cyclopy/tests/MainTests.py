import unittest
import numpy as np
import numpy.testing

class TestKernelSmoothingFunctions(unittest.TestCase):
    def setUp(self):
        self.x = np.r_[0:10:1]
        self.newx = np.r_[0:10:2]

        self.y = np.random.randn(10)
        pass

    def test_import(selfs):
        from cyclopy.NumericalMethods import KernelSmoother

    def test_create(self):
        from cyclopy.NumericalMethods import KernelSmoother
        ks = KernelSmoother(self.x, self.y)

    def test_cross_compute(self):
        from cyclopy.NumericalMethods import KernelSmoother

        ks = KernelSmoother(self.x, self.y)
        disabled = ks(self.newx, h=2, method='epa', disable_kdtree=True, just_n_count=False, just_variance=False)
        enabled = ks(self.newx, h=2, method='epa', disable_kdtree=False, just_n_count=False, just_variance=False)

        np.testing.assert_array_almost_equal(disabled, enabled)

class TestTSFunctions(unittest.TestCase):
    def setUp(self):
        self.x = np.r_[0:10:1]
        self.y = np.random.randn(10)

        from cyclopy.Elements import TimeSeriesXY, TimeSeriesEven
        self.xy = TimeSeriesXY(self.x, self.y)
        self.even = TimeSeriesEven(self.y, 1, 0.0)

    def test_resample(self):
        even = self.xy.getAsEvenSpaced(1)
        minlength = min(even.getY().size, self.even.getY().size)

        # this may fail in some cases but should be +- always true
        np.testing.assert_array_almost_equal(even.getY()[0:minlength], self.even.getY()[0:minlength], 2)

    def test_plot(self):
        self.xy.plot()
        self.even.plotMTMSpectrum()

        # self.even.plotCWT()  # disabled for now
        self.even.plotEvolutionarySpectraLines3D(self.x)
        self.even.plotHilbertTransformAbs()





if __name__ == '__main__':
    unittest.main()