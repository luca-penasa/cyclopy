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