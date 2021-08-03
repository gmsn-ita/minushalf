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

    @abstractmethod
    def get_software_name(self) -> str:
        """
        Returns the name of the software that runs first principles calculations
        """
    @abstractmethod
    def get_correction_code(self) -> str:
        """
        Returns the code used to identify the correction
        """
    @abstractmethod
    def get_overwrite_vbm(self) -> str:
        """
        Returns the parameter that overwrites vbm
        """

    @abstractmethod
    def get_overwrite_cbm(self) -> str:
        """
        Returns the parameter that overwrites cbm
        """
