"""
Parser for minushalf.yaml
"""
import yaml
from minushalf.utils import Softwares


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
        self.software = software
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

        available_softwares = Softwares.list()
        is_software_avalilable = any(element == name
                                     for element in available_softwares)
        if not is_software_avalilable:
            raise ValueError("Parameter software is not filled correctly")

        self._software = name.upper()

    @staticmethod
    def from_file(filename: str = "minushalf.yaml"):
        """
        Receives a file and catch all the parameters
        presents in the documentation
        """
        with open(filename, "r") as file:
            parsed_input = yaml.load(file, Loader=yaml.FullLoader)

        if "software" in parsed_input:
            software = parsed_input["software"]
        else:
            software = "VASP"

        if software.lower() in parsed_input:
            software_configurations = parsed_input[software.lower()]
        else:
            software_configurations = None

        if "atomic_program" in parsed_input:
            atomic_program = parsed_input["atomic_program"]
        else:
            atomic_program = None

        if "correction" in parsed_input:
            correction = parsed_input["correction"]
        else:
            correction = None

        return MinushalfYaml(
            software,
            software_configurations,
            atomic_program,
            correction,
        )
