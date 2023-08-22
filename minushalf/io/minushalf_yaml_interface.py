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

    @abstractmethod
    def get_vbm_characters(self) -> list:
        """
        Returns the parameter that cbm_characters
        """

    @abstractmethod
    def get_cbm_characters(self) -> list:
        """
        Returns the parameter vbm_characters
        """

    @abstractmethod
    def get_potential_folder(self) -> str:
        """
        Returns the potential folder name
        """

    @abstractmethod
    def get_amplitude(self) -> float:
        """
        Returns the amplitude
        """

    @abstractmethod
    def get_max_iterations(self) -> int:
        """
        Returns max iterations
        """

    @abstractmethod
    def get_exchange_corr_code(self) -> str:
        """
        Returns exchange correlation code
        """

    @abstractmethod
    def get_calculation_code(self) -> str:
        """
        Returns the calculation code
        """

    @abstractmethod
    def get_valence_cut_initial_guess(self) -> str:
        """
        Returns the valence cut initial guess
        """

    @abstractmethod
    def get_conduction_cut_initial_guess(self) -> str:
        """
        Returns the conduction cut initial guess
        """

    @abstractmethod
    def get_tolerance(self) -> float:
        """
        Returns the tolerance
        """

    @abstractmethod
    def get_indirect(self) -> bool:
        """
        Returns the indirect
        """

    @abstractmethod
    def get_divide_character(self) -> list:
        """
        Returns the divide characters
        """
