"""
List the available methods to derive cut initial guess from the
nearest neighbor distance.
"""
from enum import Enum, unique


@unique
class CutInitialGuessMethods(Enum):
    """
    Enum type for the methods of derive initial guess.
    """

    three_dimensions = "3d"  # All electrons

    def __str__(self):
        return str(self.name)

    @staticmethod
    def to_list():
        """
        Generate list of available calculation codes
        """
        return list(map(lambda element: element.value, CutInitialGuessMethods))
