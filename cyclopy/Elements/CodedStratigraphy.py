# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 18:50:51 2014

@author: luca
"""
from __future__ import division

import numpy as np


from .Stratum import Stratum

class CodedStratigraphy(object):
    """
    Is an ordered collection of Stratum.py objects that build up a stratigraphic column
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

    def update_secondary_data(self):
        self._update_limits()


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

    def get_transition_matrix(self, only_counts=False):
        number = self.get_number_of_different_lithologies()
        ids = self.get_unique_ids().tolist()

        counts = np.zeros((number, number))

        for i, s in enumerate(self):
            if i == 0:
                last = s
                continue
            else:
                counts[ids.index(s.id), ids.index(last.id)] += 1
                last = s
        if only_counts is True:
            return counts
        else:
            return (counts.T/ np.sum(counts, axis=1)).T

    def get_synthetic_series(self):
        T = self.get_transition_matrix()
        n = len(self)  # same length of this series
        stratigraphy = CodedStratigraphy()

        # start code is always the most abundant
        codes = self.get_unique_ids()

        numbers = []
        for id in codes:
            numbers.append(len(self.get_strata_of_id(id)))

        start_id = codes[np.argmax(numbers)]
        # print("using as starting id: " + str(start_id)) # was for debug only

        s = Stratum(0.0, start_id)
        stratigraphy.strata.append(s)

        old = s
        for i in range(n-1):
            # we determine which lithocode is appropiate
            prev_arg = codes.tolist().index(old.id)
            next_arg = self._do_MC_prediction(T, prev_arg)

            next_id = codes[next_arg]

            s = Stratum(0.0, next_id)
            stratigraphy.strata.append(s)

            old = s


        # now in a second pass we decide widths
        mus = []
        stds = []
        mins = []
        maxs = []
        for id in codes:
            mu, std = self.get_width_stats(id)
            min, max = self.get_width_min_max(id)
            mus.append(mu)
            stds.append(std)
            mins.append(min)
            maxs.append(max)

        for i, s in enumerate(stratigraphy):
            id = s.id
            arg_id = codes.tolist().index(id)
            thick = CodedStratigraphy._generate_gaussian_bounded(mus[arg_id], stds[arg_id], mins[arg_id], maxs[arg_id])
            stratigraphy.strata[i].thickness = thick



        stratigraphy.update_secondary_data()

        return stratigraphy

    def _do_MC_prediction(self, T, previous):
        """
        :param T: is the nxn prob matrix with the transitions probabilities
        :param previous: is the id in range(n) for the previous state
        :return: the next state as integer between 0 and n
        """
        n = T.shape[0] # n of rows that is == n cols
        selected = np.random.choice(n, p = T[previous])
        return selected

    def get_width_stats(self, id):

        wds = self.get_list_of_widths(id)

        return np.mean(wds), np.std(wds)

    def get_width_min_max(self, id):
        wds = self.get_list_of_widths(id)
        return np.min(wds), np.max(wds)

    @staticmethod
    def _generate_gaussian_bounded(mu, std, min, max):
        import random
        number = random.gauss(mu, std)

        if number < min:
            return min
        elif number > max:
            return max
        else:
            return number

    def change_litho_id(self, old_id, new_id):
        for s in self:
            if s.id ==old_id:
                s.id = new_id



    def get_shuffled_series(self, method="advanced", moveable=[]):
        """
        Get a CodedStratigraphy element with the same Strata but shuffled in some way.
        This can be used to evaluate noise level wehn performing spectral estimation on coded
        stratigraphy.
        Please remmber that cyclostratigraphy based on coded lithologs is not a really good idea in general.
        But in some cases it is the only thing you can do.

        :param method='base' or 'advanced'
            base:       -> simply shuffle the stratigraphic elements using an uniform distribution.
            advanced:   -> shuffle elements by type. keeping relationships between elemets fixed. This means we simply shuffle
                           widths of strata
            widths;     -> shuffle the widths ot strata, but only within the same strata type - should be equal to advanced,
                           but with a different implementation
            move:       -> we first use widths algorithm to produce a series with scrambled widths then we move the
                           lithologies with a code in moveable to random positions in the section.
            mc1         -> use irst order markov chain to simulate the sequence the use gaussian distribuions to chose
                           widths of strata - THIS is not shuffling.
        :param moveable a list containing the id codes of lithologies to be moved
        :raise ValueError if a wrong method is passed

        """


        from copy import deepcopy
        cop = deepcopy(self)

        if method == "simple":
            np.random.shuffle(cop.strata)

        elif method == "advanced":
            ids = self.get_unique_ids()
            for id in ids:
                posistions = self.get_strata_of_id(id, only_positions=True )
                posistions_shuffled = posistions[:]
                np.random.shuffle(posistions_shuffled)  # perform shuffling
                for orig, shuffled in zip(posistions, posistions_shuffled):
                    a_stratum = cop[orig]
                    b_stratum = cop[shuffled]

                    # swap them
                    cop[orig] = b_stratum
                    cop[shuffled] = a_stratum

        elif method == 'widths':
            for id in cop.get_unique_ids():
                thicknesses = cop.get_list_of_widths(id)
                np.random.shuffle(thicknesses)
                for s, newt in zip(cop.get_strata_of_id(id), thicknesses):
                    s.thickness = newt

        elif method == 'move':
            cop = cop.get_shuffled_series(method='widths')
            for id in moveable:
                orig_pos = cop.get_strata_of_id(id, only_positions=True)
                print(orig_pos)

                possible_positions = np.arange(len(cop))

                np.random.shuffle(possible_positions)
                print(possible_positions[0:len(orig_pos)])

                for (old_pos, new_pos) in zip(orig_pos, possible_positions[0:len(orig_pos)]):
                    print("moved from " + str(old_pos) + " to " + str(new_pos))
                    cop.strata.insert(new_pos, cop.strata.pop(old_pos))
        elif method == 'mc1':
            return self.get_synthetic_series()

        else:
            raise ValueError("wrong argument selected")


        cop.update_secondary_data()

        return cop

    def get_list_of_widths(self, id = None):
        if id == None:
            return [a.thickness for a in self ] #  return all
        else:
            return [a.thickness for a in self if (a.id == id)] # only for matches


    def get_strata_of_id(self, id, only_positions=False):
        """
        id is a stratum.id value to use to filter out strata
        it will return an empty list if those indices does not exist
        if only_positions == True it will retun the position of the requested strata in the self.strata list
        """
        if only_positions == False:
            result = []
            for s in self:
                if s.id == id:
                    result.append(s)

            return result
        else:
            return [arg for (arg, item) in enumerate(self.strata) if (item.id == id)]

    def __getitem__(self, item):
        return self.strata.__getitem__(item)

    def __setitem__(self, key, value):
        assert (isinstance(value, Stratum))
        self.strata.__setitem__(key, value)


    def get_unique_ids(self):
        """
        get the list of existent ids in the current stratigraphic log

        list will be void if th stratigraphic log is empty without any warning
        """

        return np.unique([s.id for s in self])

    def get_number_of_different_lithologies(self):
        return len(self.get_unique_ids())



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
        u_limits = self.get_upper_limits()

        ids = []
        for i, (low, up) in enumerate(zip(l_limits, u_limits)):
            if (position >= low) & (position < up):
                ids.append(i)

        if len(ids) == 0:
            print("Pertinent stratum not found")
            return None

        else:
            return self.strata[ids[0]]

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


    def to_xy_series(self, step=0.01):
        """
        Returns a time series evenly spaced obtained sampling this stratigrahic column.
        is the inverse operation of from_xy_series

        it will return two numpy arrays, x and y
        """
        x = np.r_[np.min(self.limits)+0.5*step:np.max(self.limits) - 0.5*step:step]

        y = []
        for pos in x:
            y.append(self.get_stratum_at_position(pos).id)

        return x, np.array(y)

    def to_evenly_spaced_series(self, *args, **kwargs):

        from ..Elements import TimeSeriesEven
        x, y = self.to_xy_series(*args, **kwargs)

        ser = TimeSeriesEven(y, x[1] - x[0], np.min(x))

        return ser


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




