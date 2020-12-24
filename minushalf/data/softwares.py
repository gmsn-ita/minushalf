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
    def list():
        return list(map(lambda element: element.value, Softwares))
