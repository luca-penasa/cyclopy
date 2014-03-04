__author__ = 'luca'


import numpy as np
from scipy.stats import kde

class NoiseModelerCodedStratigraphy():
    """
    Monte Carlo Noise modeling for coded stratigraphy.
    It expects as input the stratigraphy as color coded evenly-spaced time-series
    e.g. [0,0,0,1,1,1,0,0,0,2,2,2,2,2,2,1,1,0,0] is a valid array.
    """

    def __init__(self, coded_strat):
        """
        coded_strat is an array-like object of color-coded, evenly-spaced stratigraphy.
        Integers are expected (coded_strat will be forced to dtype=int)
        """
        self.input_series = np.asarray(coded_strat, dtype=int)

        # now extract data and stats from this series
        self.unique_codes = np.unique(self.input_series)

        #split up the series
        self.splits = self._split_by_codes()

        # now for each unique code count the size of each split
        self.splits_sizes = self._count_split_size_per_unique_code()

        self.splits_n = self._count_n_of_splits_per_unique_code()

        # this will count how many samples for litho-type are present in the series
        self.unique_sizes = self._compute_overall_unique_sizes()

        # establish what is the code to be used as background
        # other codes will be placed randomly in the series on top of this background
        self.background_code = np.max(self.unique_codes[self.unique_sizes.argmax()])

        #now fit a distribution for each unique code
        self.distributions = self._fit_ditributions_to_uniques()


    def _split_by_codes(self):
        """
        :rtype : is a list of numpy int-typed arrays
        """
        inds = np.where(np.diff(self.input_series))[0]
        splits = np.split(self.input_series, inds + 1)
        return splits

    def _compute_overall_unique_sizes(self):
        out = []
        for u in self.splits_sizes:
            out.append(np.sum(u))

        return np.array(out)


    def _count_split_size_per_unique_code(self):
        """
        :return: a list of sizes for each unique type in input series
        """
        codes = []
        count = []
        for s in self.splits:
            codes.append(s[0])
            count.append(len(s))

        out = []
        for code in self.unique_codes:
            ids = np.where(codes == code)[0]
            out.append(np.array(count)[ids])

        return out


    def _count_n_of_splits_per_unique_code(self):
        n = []
        for s in self.splits_sizes:
            n.append(len(s))

        return n


    def _fit_ditributions_to_uniques(self):
        """

        :rtype : a list
        """
        distributions = []
        for values in self.splits_sizes:
            distributions.append(kde.gaussian_kde(values))

        return distributions

    @staticmethod
    def _sample_distribution(distribution, sample_number, range=(None, None)):
        """
        :param distribution:
        :param sample_number:
        :param range: is a tuple defining a lower and upper bound for samples that we consider valid
                      it may be sorted or unsorted as you like (will be sorted)
        """

        # ensure range is sorted
        range = np.sort(range)

        if range[0] is None and range[1] is None:
            return distribution.resample(sample_number)[0].astype(int)

        else:
            # we sample two times the requested number
            vals = distribution.resample(sample_number * 2)[0].astype(int)
            if range[0] is not None:
                # cut out values lower than range[0]
                vals = vals[vals > range[0]]

            if range[1] is not None:
                # cut out values greater than range[1]
                vals = vals[vals < range[1]]

            # now check if we have enough values
            if len(vals) < sample_number:
                print("NOT ENOUGH SAMPLES!!!")
                return None

            else:
                return vals[0:sample_number]

        pass

    def generateRandomSeries(self):
        N = len(self.input_series)  # requested lenght
        back_code = self.background_code

        oserie = np.ones(N) * back_code

        good_uniques_ids = np.argwhere(self.unique_codes != back_code)

        for id in good_uniques_ids:
            this_code = self.unique_codes[id]  # we need to insert this id into out series
            n_elements = self.splits_n[id]  # we need to insert n_elements elements :-)
            dist = self.distributions[id]

            #get random thicknesses
            thicknesses = self._sample_distribution(dist,
                                                    n_elements,
                                                    (0, None))  #bigger than zero, but with no upper limit

            # now we insert these thicknesses at random positions, from an uniform distribution