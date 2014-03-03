__author__ = 'luca'
__doc__ = """
This module contains generic code for goeology related operations
"""

# only importing most used methods
from cyclopy.GeologyHelpers.GeometricGeology import DipDirectionAndDipAngleToDipVector, DipDirectionToStrikeVector, \
    Project3DPointsOnStratigraphicNormal

from cyclopy.GeologyHelpers.StratigraphicModeling import OutcropModel, UndeformedModel, ResidualsModel