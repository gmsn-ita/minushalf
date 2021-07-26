"""
Class for atomic program software parameters in minushalf.yaml
"""
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
