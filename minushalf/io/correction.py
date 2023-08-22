"""
Class for correction input parameters in minushalf.yaml
"""
from shutil import Error
from minushalf.io.minushalf_yaml_tags_interface import MinushalfYamlTags
import loguru
from minushalf.io.correction_code import CorrectionCode
from minushalf.io.minushalf_yaml_tags_interface import MinushalfYamlTags


class Correction(MinushalfYamlTags):
    """
    Set parameters and their default values
    """

    def __init__(self,
                 correction_code: str = CorrectionCode.get_default(),
                 potfiles_folder: str = "minushalf_potfiles",
                 amplitude: float = 1.0,
                 valence_cut_guess: list = None,
                 conduction_cut_guess: list = None,
                 tolerance: float = 0.01,
                 fractional_valence_treshold: float = 10,
                 fractional_conduction_treshold: float = 9,
                 cbm_characters: list = None,
                 vbm_characters: list = None,
                 overwrite_vbm: list = None,
                 overwrite_cbm: list = None,
                 indirect: bool = False,
                 divide_character: list = None) -> None:
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
                indirect: Realize calcualtions indirect
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
        self.cbm_characters = cbm_characters
        self.vbm_characters = vbm_characters
        self.indirect = indirect
        self.divide_character = divide_character

    @property
    def divide_character(self) -> list:
        """
        Returns:
            Factor that divides the correction between atoms
        """
        return self._divide_character

    @divide_character.setter
    def divide_character(self, factors: int) -> None:
        """
        Set factor that divides the correction between atoms
        """
        try:
            if factors != None:
                for element in factors:
                    element[0] = element[0].capitalize()
                    element[1] = element[1].lower()
                    element[2] = int(element[2])
        except:
            loguru.logger.error("divide_character incorrectly specified")
            raise Error("divide_character incorrectly specified")

        self._divide_character = factors

    @property
    def vbm_characters(self) -> list:
        """
        Returns:
            Artificial character in vbm
        """
        return self._vbm_characters

    @vbm_characters.setter
    def vbm_characters(self, characters: list) -> None:
        """
        Set vbm_characters 
        """
        try:
            if characters != None:
                for element in characters:
                    element[0] = element[0].capitalize()
                    element[1] = element[1].lower()
                    element[2] = float(element[2])
        except:
            loguru.logger.error("repalce_vbm incorrectly specified")
            raise Error("repalce_vbm incorrectly specified")

        self._vbm_characters = characters

    @property
    def cbm_characters(self) -> list:
        """
        Returns:
            Artificial character in cbm
        """
        return self._cbm_characters

    @cbm_characters.setter
    def cbm_characters(self, characters: list) -> None:
        """
        Set cbm_characters 
        """

        try:
            if characters != None:
                for element in characters:
                    element[0] = element[0].capitalize()
                    element[1] = element[1].lower()
                    element[2] = float(element[2])
        except:
            loguru.logger.error("repalce_cbm incorrectly specified")
            raise Error("repalce_cbm incorrectly specified")

        self._cbm_characters = characters

    @property
    def valence_cut_guess(self) -> list:
        """
        Returns:
            CUT guess for nelder mead algorithm
        """
        return self._valence_cut_guess

    @valence_cut_guess.setter
    def valence_cut_guess(self, cut_guess: list) -> None:
        """
        Set valence_cut_guess
        """

        try:
            if cut_guess != None:
                for element in cut_guess:
                    element[0] = element[0].capitalize()
                    element[1] = element[1].lower()
                    element[2] = float(element[2])
        except:
            loguru.logger.error("valence_cut_guess incorrectly specified")
            raise Error("valence_cut_guess incorrectly specified")

        self._valence_cut_guess = cut_guess

    @property
    def conduction_cut_guess(self) -> list:
        """
        Returns:
            CUT guess for nelder mead algorithm
        """
        return self._conduction_cut_guess

    @conduction_cut_guess.setter
    def conduction_cut_guess(self, cut_guess: list) -> None:
        """
        Set conduction_cut_guess
        """

        try:
            if cut_guess != None:
                for element in cut_guess:
                    element[0] = element[0].capitalize()
                    element[1] = element[1].lower()
                    element[2] = float(element[2])
        except:
            loguru.logger.error("conduction_cut_guess incorrectly specified")
            raise Error("conduction_cut_guess incorrectly specified")

        self._conduction_cut_guess = cut_guess

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
        parameters_dict = self.__dict__.copy()

        # removing private variables
        parameters_dict["correction_code"] = parameters_dict.pop(
            "_correction_code", None)
        parameters_dict["overwrite_vbm"] = parameters_dict.pop(
            "_overwrite_vbm", None)
        parameters_dict["overwrite_cbm"] = parameters_dict.pop(
            "_overwrite_cbm", None)
        parameters_dict["vbm_characters"] = parameters_dict.pop(
            "_vbm_characters", None)
        parameters_dict["cbm_characters"] = parameters_dict.pop(
            "_cbm_characters", None)
        parameters_dict["divide_character"] = parameters_dict.pop(
            "_divide_character", None)
        parameters_dict["valence_cut_guess"] = parameters_dict.pop(
            "_valence_cut_guess", None)
        parameters_dict["conduction_cut_guess"] = parameters_dict.pop(
            "_conduction_cut_guess", None)

        return parameters_dict
