"""
Interface for minushalf.yaml
"""
from abc import ABC, abstractmethod


class MinushalfYaml(ABC):
    """
    Interface
    """
    @abstractmethod
    def get_correction_params(self) -> dict:
        """
        Get dictionary of correction parameters
        """

    @abstractmethod
    def get_atomic_program_params(self) -> dict:
        """
        Get dictionary of atomic program parameters
        """

    @abstractmethod
    def get_software_configurations_params(self) -> dict:
        """
        Get dictionary of software configurations parameters
        """
