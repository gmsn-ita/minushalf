"""
List calculation code options for the INP file
"""
from enum import Enum, unique


@unique
class CalculationCode(Enum):
    """
    Enum type for the calculation code
    """

    ae = "ae"  # All electrons

    def __str__(self):
        return str(self.name)

    @staticmethod
    def to_list():
        """
        Generate list of available calculation codes
        """
        return list(map(lambda element: element.value, CalculationCode))
