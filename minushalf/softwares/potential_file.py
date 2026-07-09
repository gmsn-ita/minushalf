"""
Interface for
classes that handles with the 
local potential
"""
from abc import ABC, abstractmethod


class PotentialFile(ABC):
    """
    Interface for defining common methods
    between classes that handles with potential
    files of different softwares
    """
    @abstractmethod
    def to_stringlist(self) -> list:
        """
        Abstract methods that returns a list
        containing the lines of the potential file
        """

    @abstractmethod
    def to_file(self, filename: str) -> None:
        """
        Abstract methods that output
        a potential file with a specific
        filename
        """

    @abstractmethod
    def get_local_potential(self) -> list:
        """
        Abstract methods returns
        the local potential
        """

    @abstractmethod
    def get_maximum_module_wave_vector(self) -> float:
        """
        Abstract methods returns
        the maximum modulus of
        the wave vector in reciprocal space
        """

    @abstractmethod
    def get_name(self) -> float:
        """
        Abstract methods returns
        the name of the potential file
        """
