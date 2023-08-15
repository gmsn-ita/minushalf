"""
List softwares supported
by the CLI
"""
from enum import Enum, unique
from minushalf.softwares.vasp_factory import Vasp


@unique
class Softwares(Enum):
    """
    Enum type for the softwares supported by the program
    """

    vasp = "VASP"

    def __str__(self):
        return str(self.name)

    @staticmethod
    def get_default():
        """
        Returns the default value for this parameter
        """
        return Softwares.vasp.value

    @staticmethod
    def to_list():
        """
        Generate list of available softwares
        """
        return list(map(lambda element: element.value, Softwares))


def get_software_factory(software: str):
    software_to_factory = {
        Softwares.vasp.value: Vasp()
    }
    return software_to_factory[software]
