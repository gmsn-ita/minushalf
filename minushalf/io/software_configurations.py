"""
Class for atomic program software parameters in minushalf.yaml
"""
import loguru
from minushalf.data import Softwares
from minushalf.interfaces import MinushalfYamlTags


class SoftwareConfigurations(MinushalfYamlTags):
    """
    Set parameters and their default values
    """
    def __init__(
            self,
            command: list = None,
            software_name: str = Softwares.get_default(),
    ) -> None:
        """
            Args:
                command (list):  Command to perform first-principles calculations
                software_name (str): Name of software that performs first principles calculations
        """
        self.software_name = software_name
        self.command = command

    @property
    def command(self) -> dict:
        """
        Returns:
            Command to perform first-principles calculations
        """
        return self._command

    @command.setter
    def command(self, command_list: list) -> None:
        """
        Set command variable
        """
        if not command_list:
            software_commands = {Softwares.vasp.value: ['mpirun', 'vasp']}
            command_list = software_commands[self.software_name]

        self._command = command_list

    @property
    def software_name(self) -> str:
        """
        Returns:
            Name of the software used for ab initio calculations (VASP,...)
        """
        return self._software_name

    @software_name.setter
    def software_name(self, name: str) -> None:
        """
        Verify if the software name is a valid name

        Args:
            name (str): Name of the software
        """

        available_softwares = Softwares.to_list()
        is_software_avalilable = any(element.lower() == name.lower()
                                     for element in available_softwares)
        if not is_software_avalilable:
            loguru.logger.error("Parameter software is not filled correctly")
            raise ValueError("Parameter software is not filled correctly")

        self._software_name = name.upper()

    def to_list(self):
        """
        return list with the class variables
        """
        return [self._command]

    def to_dict(self):
        """
        Return dictionary with the class variables
        """
        return {"command": self._command}
