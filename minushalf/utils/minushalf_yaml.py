"""
Parser for minushalf.yaml
"""
import yaml
from loguru import logger
from minushalf.data import (
    Softwares,
    VaspDefaultParams,
    AtomicProgramDefaultParams,
    CorrectionDefaultParams,
    CorrectionCode,
    MinushalfParams,
)


class MinushalfYaml():
    """
    Class that parses the input
    for the execute command
    """
    def __init__(
        self,
        software: str,
        software_configurations: dict,
        atomic_program: dict,
        correction: dict,
    ):
        """
        Constructs a class for input in the execute command
        """
        self.software = software.upper()
        self.software_configurations = software_configurations
        self.atomic_program = atomic_program
        self.correction = correction

    @property
    def software(self) -> str:
        """
        Returns:
            Name of the software used for ab initio calculations (VASP,...)
        """
        return self._software

    @software.setter
    def software(self, name: str) -> None:
        """
        Verify if the symbol is a valid periodic table element and
        format the string correctly.

        Args:
            symbol (str): chemical symbol of the element (H, He, Li...)
        """

        available_softwares = Softwares.to_list()
        is_software_avalilable = any(element == name
                                     for element in available_softwares)
        if not is_software_avalilable:
            logger.error("Parameter software is not filled correctly")
            raise ValueError("Parameter software is not filled correctly")

        self._software = name.upper()

    @property
    def software_configurations(self) -> str:
        """
        Returns:
            Some configurations for the softwares that runs
            ab initio calculations
        """
        return self._software_configurations

    @software_configurations.setter
    def software_configurations(self, configurations: dict) -> None:
        """
            Verify if the file has all the cofigurations needed. If not, use
            the default values
                Args: configurations (dict): Configuration parameters for each software
                    VASP:
                        number_of_cores: Number of cores to run the proccess. The
                        default is one.
                        path_to_vasp: Path to the VASP software.The deafault is 'vasp'
        """
        choice_params = {Softwares.vasp.value: VaspDefaultParams}

        software_params = choice_params[self.software]
        default_params = software_params.to_dict()

        if not configurations:
            configurations = default_params
        else:
            for name, value in default_params.items():
                if not name in configurations:
                    configurations[name] = value
            for name, value in configurations.items():
                if not name in default_params:
                    raise ValueError("param {} is not defined".format(name))

        self._software_configurations = configurations

    @property
    def atomic_program(self) -> dict:
        """
        Returns:
            Some configurations for the atomic program
        """
        return self._atomic_program

    @atomic_program.setter
    def atomic_program(self, configurations: dict) -> None:
        """
        Verify if the file has all the cofigurations needed. If not, use
        the default values

            Args:
                configurations(dict): Configurations params for
                atomc program. Contain the following fields:

                    exchange_correlation_type: Default is 'pb'

                    calculation_code: Default is 'ae'

                    max_iterations: Default is 100
        """
        default_params = AtomicProgramDefaultParams.to_dict()
        if not configurations:
            configurations = default_params
        else:
            for name, value in default_params.items():
                if not name in configurations:
                    configurations[name] = value
            for name, value in configurations.items():
                if not name in default_params:
                    raise ValueError("param {} is not defined".format(name))

        self._atomic_program = configurations

    @property
    def correction(self) -> dict:
        """
        Returns:
            Some configurations for the correction in the
            potential files
        """
        return self._correction

    @correction.setter
    def correction(self, configurations: dict) -> None:
        """
        Verify if the file has all the cofigurations needed. If not, use
        the default values

            Args:
                configurations(dict): Configurations params for
                correct potential. Contain the following fields:

                    correction_code: Default is 'v'

        """
        default_params = CorrectionDefaultParams.to_dict()
        if not configurations:
            configurations = default_params
        else:
            for name, value in default_params.items():
                if not name in configurations:
                    configurations[name] = value
            for name, value in configurations.items():
                if not name in default_params:
                    raise ValueError("param {} is not defined".format(name))

        available_correction_codes = CorrectionCode.to_list()
        is_code_available = any(element == configurations[
            CorrectionDefaultParams.correction_code.__str__()]
                                for element in available_correction_codes)
        if not is_code_available:
            raise ValueError("Invalid value for correction code")

        self._correction = configurations

    @staticmethod
    def from_file(filename: str = "minushalf.yaml"):
        """
        Receives a file and catch all the parameters
        presents in the documentation
        """
        parsed_input = {}
        try:
            with open(filename, "r") as file:
                parsed_input = yaml.load(file, Loader=yaml.FullLoader)
        except FileNotFoundError:
            logger.info("File not found, default parameters will be used")

        if MinushalfParams.software.value in parsed_input:
            software = parsed_input[MinushalfParams.software.value]
        else:
            software = Softwares.vasp.value

        if software.lower() in parsed_input:
            software_configurations = parsed_input[software.lower()]
        else:
            software_configurations = None

        if MinushalfParams.atomic_program.value in parsed_input:
            atomic_program = parsed_input[MinushalfParams.atomic_program.value]
        else:
            atomic_program = None

        if MinushalfParams.correction.value in parsed_input:
            correction = parsed_input[MinushalfParams.correction.value]
        else:
            correction = None

        return MinushalfYaml(
            software,
            software_configurations,
            atomic_program,
            correction,
        )
