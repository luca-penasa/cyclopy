# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 18:50:51 2014

@author: luca
"""
from __future__ import division

import numpy as np

from scipy.stats import kde



class Stratum(object):
    """
    Is a object that representa geoogical stratum, which type is identified by and id in lit_code
    """
    def __init__(self, thickness=None, lit_code=None):
        self.id = int(lit_code)  # only integers codes are right
        self.thickness = float(thickness)

    @property
    def thickness(self):
        return self._thickness

    @thickness.setter
    def thickness(self, value):
        self._thickness = value


    @property
    def id(self):
        return self._id

    @id.setter
    def id(self, value):
        self._id = value


    def __repr__(self):
        return str(self.id) + " T: " + str(self.thickness)

    def sum_thickness(self, other):
        return self.thickness + other.thickness


class CodedStratigraphy(object):
    """
    Is an ordered collection of Stratum objects that build up a stratigraphic column
    Altohugh it may not be the more effective way to represent a stratigraphic column
    it is a well accepted way in geology. The first stratum will have a positioning in the
    stratigraphic domain (a stratigraphic position or SP) that is the base of this stratum.
    Strata are expected to be ordered from lower SP to higher SP.


    methods for inserting stratum are:
    insert_new_stratum      - insert a stratum in given position, resize others
    push_strata             - push a new strata to the db


    properties:
    limits                  - the upper and lower limits for each strata
    """
    def __init__(self):
        self._strata = list()
        self._lower_position = 0.0

        self._limits = np.array([])

    @property
    def limits(self):
        return self._limits


    def __repr__(self):
        row_format ="ID: {:<4d}, Thickness: {:<6.2f}  \n"
        string = str()
        for s in self.strata[::-1]:
            string += row_format.format(s.id,s.thickness)

        return string


    @property
    def strata(self):
        return self._strata

    @property
    def lower_position(self):
        return self._lower_position

    @lower_position.setter
    def lower_position(self, val):
        #todo add type checking please
        self._lower_position = val
        self._update_limits()


    def push_strata(self, strata):
        assert isinstance(strata, Stratum)
        self.strata.append(strata)
        self._update_limits()

    def _update_limits(self):
        self._limits = np.zeros(len(self) + 1)

        current = self.lower_position

        for i,s in enumerate(self):
            self._limits[i] = current
            current += s.thickness

        self._limits[-1] = current


    def get_limits_for_stratum(self, stratum):
        """
        Return lower and upper limits for the strata with the given ID
        """
        assert isinstance(stratum, Stratum)

        try:
            id = self.strata.index(stratum)

        except ValueError:
            raise ValueError("This strata do not belong to this stratigraphic log")

        return self.limits[id], self.limits[id + 1]


    def merge_consecutive_equal_strata(self):
        out = []

        for i, s in enumerate(self.strata):
            if i == 0:
                cur_strat = s
                continue

            if s.id  == cur_strat.id:  # are of the same type
                cur_strat.thickness += s.thickness

            else:  # they are not of the same type
                out.append(cur_strat)
                cur_strat = s

            if i == len(self.strata) - 1:
                    out.append(cur_strat)

        self.strata[:] = out[:]



    def __iter__(self):
        return iter(self.strata)

    def __len__(self):
        return len(self.strata)



    @staticmethod
    def from_xy_series(x,y):
        assert(len(x) == len(y))

        sort_ids = np.argsort(x)
        x = x[sort_ids]
        y = y[sort_ids]

        assert(len(np.unique(np.round(np.diff(x), 4))) == 1)  # we need to round them

        step = np.diff(x)[0]

        column = CodedStratigraphy()

        for i in y:
            column.push_strata(Stratum(step, i))

        column.lower_position = np.min(x)

        column.merge_consecutive_equal_strata()
        column._update_limits()

        return column

    def is_stratum_in_range(self, stratum, range):
        """
        It will tell you if a given strata fall within the given range of stratigraphic positions
        aka given a range of stratigraphic positions it will say if your strata is within it,
        if it is only partially it will gives True
        """

        assert(isinstance(stratum, Stratum))

        low, up = self.get_limits_for_stratum(stratum)


    def get_stratum_at_position(self, position):
        l_limits = self.get_lower_limits()
        id = [i for i in l_limits if i > position][0]
        return self.strata[id]

    def get_slice(self, lowSP, highSP):
        """
        Return a subset of the stratiraphy with the interested strata
        Partially comprised strata are included
        :param lowSP:
        :param highSP:
        """
        assert(lowSP < highSP)

        # now we use the limits to locate interested strata
        # good_ids = \

        args = np.argwhere ( (self.get_upper_limits() < highSP) & (self.get_lower_limits() > lowSP) ).T[0]



        print(self.get_upper_limits() < highSP)

        print(self.get_lower_limits() > lowSP)
        print (args)

        # return [s for s in self.strata if ]

    def get_upper_limits(self):
        return self.limits[1::]

    def get_lower_limits(self):
        return self.limits[:-1]




    def insert_new_stratum(self, stratum, position, method='smart'):
        """
        inserts a new sratum in the position given by position
        the position is intended o be the lower position of the straturm. (not the center or upper)
        method controls the insertion policies
        methods are:
        - 'smart' it will merge the given strata ONLY if it will not completely erease any other yet-existent strata

        returns:
            True if insertion was sudccesful
            False if it the stratum cannot be inserted with the give method
        """
        interested_strata = self.get_slice(position, position + stratum.thickness)

        print(interested_strata)

        if method == 'smart':
            if len(interested_strata) >= 3:
                print("Cannot perform merging. It would erease 1 or more strata.")
                return False

            else:
                # we resize the two interested elements
                s1 = interested_strata[0]
                s2 = interested_strata[1]

                L1, L2 = self.get_limits_for_stratum(s1) # upper and lower limits for the lower interested stratum
                L3 = position # lower limit for the new strata to insert


                T1 = interested_strata[0].thickness  # thickness of lower interested stratum
                T2 = interested_strata[1].thickness  # thickness of upper interested stratum
                T3 = stratum.thickness  # the thicknss of the newly inserted stratum

                newT1 = T1 - (L2 - L3)  # the new thickness of the lower startum
                newT2 = (L2 + T2) - (L3 + T3)  # new thickness for the upper stratum

                # adjust the thicknesses
                s1.thickness = newT1
                s2.thichness = newT2

                # and insert the new stratum
                self.strata.insert(stratum, self.strata.index(s1) )

                # and recall an udpate for the limits
                self._update_limits()









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











