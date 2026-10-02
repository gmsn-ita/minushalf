"""
Class for atomic program software parameters in minushalf.yaml
"""
import loguru
from minushalf.softwares.softwares import Softwares
from minushalf.io.minushalf_yaml_tags_interface import MinushalfYamlTags


class SoftwareConfigurations(MinushalfYamlTags):
    """
    Set parameters and their default values.

    VASP uses a single command to perform the ab initio calculation.

    Quantum ESPRESSO (QE) uses separate commands for each executable
    involved in the workflow:
        - pw.x
        - projwfc.x
        - scf input file
    """

    def __init__(
            self,
            command: list = None,
            software_name: str = Softwares.get_default(),
            pw_command: list = None,
            projwfc_command: list = None,
            ld1_command: list = None,
            virtual_v2_command: list = None,
            input_file: str = None,
    ) -> None:
        """
        Args:
            command (list):
                VASP command used to perform first-principles calculations.

            software_name (str):
                Name of software that performs first-principles calculations.

            pw_command (list):
                QE pw.x command.

            projwfc_command (list):
                QE projwfc.x command.

            input_file (str):
                Input files needed for calculation (e.g. INCAR for VASP, input of pw.x for QE)
        """
        self.software_name = software_name

        # VASP-specific parameter
        self.command = command

        # QE-specific parameters
        self.pw_command = pw_command
        self.projwfc_command = projwfc_command
        self.ld1_command = ld1_command
        self.virtual_v2_command = virtual_v2_command

        # shared parameter
        self.input_file = input_file

    @property
    def command(self) -> list:
        """
        VASP-specific command.

        Returns:
            Command to perform first-principles calculations.
        """
        return self._command

    @command.setter
    def command(self, command_list: list) -> None:
        """
        Set VASP command.

        If no command is provided, use the default VASP command.
        """
        if not command_list:
            software_commands = {
                Softwares.vasp.value: ['mpirun', 'vasp']
            }

            # Only VASP has a default `command`.
            if self.software_name == Softwares.vasp.value:
                command_list = software_commands[self.software_name]

        self._command = command_list

    @property
    def software_name(self) -> str:
        """
        Returns:
            Name of the software used for ab initio calculations.
        """
        return self._software_name

    @software_name.setter
    def software_name(self, name: str) -> None:
        """
        Verify if the software name is valid.

        Args:
            name (str): Name of the software.
        """

        available_softwares = Softwares.to_list()
        is_software_available = any(
            element.lower() == name.lower()
            for element in available_softwares
        )

        if not is_software_available:
            loguru.logger.error(
                "Parameter software is not filled correctly"
            )
            raise ValueError(
                "Parameter software is not filled correctly"
            )

        self._software_name = name.upper()

    def to_list(self):
        """
        Return list with the class variables.
        """

        if self.software_name == Softwares.vasp.value:
            # VASP-specific
            return [self._command,
                    self.input_file]

        elif self.software_name == Softwares.qe.value:
            # QE-specific
            return [
                self.pw_command,
                self.projwfc_command,
                self.ld1_command,
                self.virtual_v2_command,
                self.input_file,
            ]

    def to_dict(self):
        """
        Return dictionary with the class variables.

        The dictionary keys must match the arguments expected by
        the corresponding software factory.
        """

        if self.software_name == Softwares.vasp.value:
            # VASP-specific
            return {
                "command": self._command,
                "input_file": self.input_file,
            }

        elif self.software_name == Softwares.qe.value:
            # QE-specific
            return {
                "pw_command": self.pw_command,
                "projwfc_command": self.projwfc_command,
                "ld1_command": self.ld1_command,
                "virtual_v2_command": self.virtual_v2_command,
                "input_file": self.input_file,
            }