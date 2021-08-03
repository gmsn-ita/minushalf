"""
Parser for minushalf.yaml
"""
import yaml
import collections
import loguru
from minushalf.data import Softwares
from .software_configurations import SoftwareConfigurations
from .atomic_program import AtomicProgram
from .correction import Correction
from minushalf.interfaces import MinushalfYamlTags,MinushalfYaml


class MinushalfYaml(MinushalfYaml):
    """
    Class that parses the input
    for the execute command
    """
    def __init__(
        self,
        software_configurations: MinushalfYamlTags,
        atomic_program: MinushalfYamlTags,
        correction: MinushalfYamlTags,
    ):
        """
        Constructs a class for input in the execute command
        """
        self.software_configurations = software_configurations
        self.atomic_program = atomic_program
        self.correction = correction

    def get_correction_params(self) -> dict:
        """
        Get dictionary of correction parameters
        """
        return self.correction.to_dict()

    def get_atomic_program_params(self) -> dict:
        """
        Get dictionary of atomic program parameters
        """
        return self.atomic_program.to_dict()

    def get_software_configurations_params(self) -> dict:
        """
        Get dictionary of software configurations parameters
        """
        return self.software_configurations.to_dict()

    def get_software_name(self) -> str:
        """
        Returns the name of the software that runs first principles calculations
        """
        return self.software_configurations.software_name

    def get_command(self) -> list:
        """
        Returns the command that runs first principles calculations
        """
        return self.software_configurations.command

    def get_correction_code(self) -> list:
        """
        Returns the code used to identify the correction
        """
        return self.correction.correction_code

    def get_overwrite_vbm(self) -> list:
        """
        Returns the parameter that overwrites vbm
        """
        return self.correction.overwrite_vbm

    def get_overwrite_cbm(self) -> list:
        """
        Returns the parameter that overwrites cbm
        """
        return self.correction.overwrite_cbm

    @staticmethod
    def _read_yaml(filename: str) -> collections.defaultdict:
        """
        Read yaml file and export a dictionary
        """
        yaml_file = {}
        try:
            with open(filename, "r") as file:
                yaml_file = yaml.load(file, Loader=yaml.FullLoader)
        except FileNotFoundError:
            loguru.logger.warning(
                "File not found, default parameters will be used")

        return collections.defaultdict(dict, yaml_file)

    @staticmethod
    def _get_default_software_name(yaml_file: collections.defaultdict) -> str:
        """
        Name of the software that performs first priniple calculations
        """
        default = Softwares.vasp.value
        name = yaml_file["software"]

        software_name = (default if name == {} else name)
        return software_name

    @staticmethod
    def from_file(filename: str = "minushalf.yaml"):
        """
        Receives a file and catch all the parameters
        presents in the documentation
        """
        yaml_file = MinushalfYaml._read_yaml(filename)

        ## Get constructor parameters
        try:
            software_name = MinushalfYaml._get_default_software_name(yaml_file)

            software_configurations = SoftwareConfigurations(
                software_name=software_name, **yaml_file[software_name.lower()])
            correction = Correction(**yaml_file["correction"])
            atomic_program = AtomicProgram(**yaml_file["atomic_program"])
        except:
            loguru.logger.error("Invalid parameters in minushalf.yaml file")
            raise KeyError("Invalid parameters in minushalf.yaml file")

        return MinushalfYaml(
            software_configurations,
            atomic_program,
            correction,
        )