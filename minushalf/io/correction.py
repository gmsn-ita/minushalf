"""
Class for correction input parameters in minushalf.yaml
"""
from minushalf.interfaces.minushalf_yaml_tags import MinushalfYamlTags
import loguru
from minushalf.data import CorrectionCode
from minushalf.interfaces import MinushalfYamlTags


class Correction(MinushalfYamlTags):
    """
    Set parameters and their default values
    """
    def __init__(self,
                 correction_code: str = CorrectionCode.get_default(),
                 potfiles_folder: str = "minushalf_potfiles",
                 amplitude: float = 1.0,
                 valence_cut_guess: float = None,
                 conduction_cut_guess: float = None,
                 tolerance: float = 0.01,
                 fractional_valence_treshold: float = 10,
                 fractional_conduction_treshold: float = 9,
                 overwrite_vbm: list = None,
                 overwrite_cbm: list = None,
                 inplace: bool = False) -> None:
        """
            Args:
                correction_code (str): Code for DFT -1/2 correction
                potfiles_folder (str): Potential files folder
                amplitude (float): Scaling factor in trimming function
                valence_cut_guess (float): Initial guess to valence correction
                conduction_cut_guess (float): Intitial guess to conduction correction
                tolerance (float): Absolute tolerance in the Nelder-Mead algorithm
                fractional_valence_treshold (float): Treshold for fractional valence correction
                fractional_conduction_treshold (float): Treshold for fractional conduction correction
                overwrite_vbm: Tag to overwrite vbm character
                overwrite_cbm: Tag to overwrite cbm character
                inplace: Realize calcualtions inplace
        """
        self.correction_code = correction_code
        self.potfiles_folder = potfiles_folder
        self.amplitude = amplitude
        self.valence_cut_guess = valence_cut_guess
        self.conduction_cut_guess = conduction_cut_guess
        self.tolerance = tolerance
        self.fractional_valence_treshold = fractional_valence_treshold
        self.fractional_conduction_treshold = fractional_conduction_treshold
        self.overwrite_cbm = overwrite_cbm
        self.overwrite_vbm = overwrite_vbm
        self.inplace = inplace

    @property
    def correction_code(self) -> dict:
        """
        Returns:
            Code for DFT -1/2 correction
        """
        return self._correction_code

    @correction_code.setter
    def correction_code(self, code: str) -> None:
        """
        Set correction_code variable
        """
        available_correction_codes = CorrectionCode.to_list()
        is_code_available = any(test_code == code
                                for test_code in available_correction_codes)
        if not is_code_available:
            loguru.logger.error(
                "Invalid value for correction code in minushalf.yaml")
            raise ValueError("Invalid value for correction code")

        self._correction_code = code

    @property
    def overwrite_cbm(self) -> dict:
        """
        Returns:
            Tag to overwrite vbm character
        """
        return self._overwrite_cbm

    @overwrite_cbm.setter
    def overwrite_cbm(self, band_location: list) -> None:
        """
        Set overwirte_cbm variable
        """
        if not band_location:
            band_location = []
        self._overwrite_cbm = band_location

    @property
    def overwrite_vbm(self) -> dict:
        """
        Returns:
            Tag to overwrite vbm character
        """
        return self._overwrite_vbm

    @overwrite_vbm.setter
    def overwrite_vbm(self, band_location: list) -> None:
        """
        Set overwirte_vbm variable
        """
        if not band_location:
            band_location = []
        self._overwrite_vbm = band_location

    def to_list(self):
        """
        return list with the class variables
        """
        return list(self.__dict__.values())

    def to_dict(self):
        """
        Return dictionary with the class variables
        """
        parameters_dict = self.__dict__

        ## removing private variables
        parameters_dict["correction_code"] = parameters_dict.pop(
            "_correction_code", None)
        parameters_dict["overwrite_vbm"] = parameters_dict.pop(
            "_overwrite_vbm", None)
        parameters_dict["overwrite_cbm"] = parameters_dict.pop(
            "_overwrite_cbm", None)

        return self.__dict__
