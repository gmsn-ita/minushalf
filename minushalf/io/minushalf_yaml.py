"""
Parser for minushalf.yaml
"""
import yaml
import collections
import loguru
from minushalf.softwares.softwares import Softwares
from minushalf.io.software_configurations import SoftwareConfigurations
from minushalf.io.atomic_program import AtomicProgram
from minushalf.io.correction import Correction
from minushalf.io.minushalf_yaml_interface import MinushalfYaml
from minushalf.io.minushalf_yaml_tags_interface import MinushalfYamlTags


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

    def get_vbm_characters(self) -> list:
        """
        Returns the parameter that cbm_characters
        """
        return self.correction.vbm_characters

    def get_cbm_characters(self) -> list:
        """
        Returns the parameter vbm_characters
        """
        return self.correction.cbm_characters

    def get_potential_folder(self) -> str:
        """
        Returns the potential folder name
        """
        return self.correction.potfiles_folder

    def get_amplitude(self) -> float:
        """
        Returns the amplitude
        """
        return self.correction.amplitude

    def get_max_iterations(self) -> int:
        """
        Returns max iterations
        """
        return self.atomic_program.max_iterations

    def get_exchange_corr_code(self) -> str:
        """
        Returns exchange correlation code
        """
        return self.atomic_program.exchange_correlation_code

    def get_calculation_code(self) -> str:
        """
        Returns the calculation code
        """
        return self.atomic_program.calculation_code

    def get_valence_cut_initial_guess(self) -> str:
        """
        Returns the valence cut initial guess
        """
        return self.correction.valence_cut_guess

    def get_conduction_cut_initial_guess(self) -> str:
        """
        Returns the conduction cut initial guess
        """
        return self.correction.conduction_cut_guess

    def get_tolerance(self) -> float:
        """
        Returns the tolerance
        """
        return self.correction.tolerance

    def get_inplace(self) -> bool:
        """
        Returns the inplace
        """
        return self.correction.inplace

    def get_divide_character(self) -> list:
        """
        Returns the divide characters
        """
        return self.correction.divide_character

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