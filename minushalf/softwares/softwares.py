"""
List softwares supported
by the CLI
"""
from enum import Enum, unique
from minushalf.softwares.vasp.vasp_factory import Vasp
from minushalf.softwares.qe.qe_factory import QE


@unique
class Softwares(Enum):
    """
    Enum type for the softwares supported by the program
    """

    vasp = "VASP"
    qe = "QE"

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
        Softwares.vasp.value: Vasp(),
        Softwares.qe.value: QE()
    }
    return software_to_factory[software]
