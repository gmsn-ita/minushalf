"""
List softwares suported
by the CLI
"""
from enum import Enum, unique


@unique
class Softwares(Enum):
    """
    Enum type for the softwares suport by the program
    """

    vasp = "VASP"

    def __str__(self):
        return str(self.name)

    @staticmethod
    def to_list():
        """
        Generate list of available softwares
        """
        return list(map(lambda element: element.value, Softwares))
