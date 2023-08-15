"""
List the available methods to obtain cut initial guess from the
nearest neighbor distance.
"""
from enum import Enum, unique


@unique
class CutInitialGuessMethods(Enum):
    """
    Enum for different methods of obtaining the initial cut guess.
    """

    three_dimensions = "3d"  # three dimensional crystals

    def __str__(self):
        return str(self.name)

    @staticmethod
    def to_list():
        """
        Generate list of available methods
        """
        return list(map(lambda element: element.value, CutInitialGuessMethods))
