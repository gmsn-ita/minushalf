"""
Implementation of the runner interface for Quantum ESPRESSO.

Provides the concrete runner class responsible for executing
the Quantum ESPRESSO workflow.
"""

import re
import shlex
import subprocess
from pathlib import Path
from typing import List

from minushalf.softwares.runner import Runner


class QERunner(Runner):
    """
    Runner for Quantum ESPRESSO calculations.

    Executes the QE workflow:

        1. pw.x < scf.in
        2. Create proj.in
        3. projwfc.x < proj.in
    """

    def __init__(
        self,
        pw_command: List[str],
        projwfc_command: List[str],
        input_file: str,
    ):
        """
        Args:
            pw_command: Command used to execute pw.x.
            projwfc_command: Command used to execute projwfc.x.
            input_file: Path to the SCF input file.
        """
        self.pw_command = pw_command
        self.projwfc_command = projwfc_command
        self.input_file = input_file

    def run(self, cwd: str = "."):
        """
        Run the Quantum ESPRESSO calculation workflow.
        """
        self._run_pw(cwd)
        self._create_proj_input(cwd)
        self._run_projwfc(cwd)

    def _run_pw(self, cwd: str = "."):
        """
        Run pw.x using the SCF input file.
        """

        scf_filename = Path(self.input_file).name

        command = (
            f"{shlex.join(self.pw_command)} "
            f"< {shlex.quote(scf_filename)}"
        )

        subprocess.run(
            command,
            shell=True,
            check=True,
            cwd=cwd,
        )

    def _create_proj_input(self, cwd: str = "."):
        """
        Create the proj.in file for projwfc.x.

        The `prefix` and `outdir` values are extracted from the SCF input
        file. The remaining PROJWFC parameters are fixed constants.
        """

        scf_filename = Path(self.input_file).name
        scf_path = Path(cwd) / scf_filename

        with open(scf_path, "r") as scf_file:
            scf_content = scf_file.read()

        prefix_match = re.search(
            r"\bprefix\s*=\s*['\"]([^'\"]+)['\"]",
            scf_content,
            re.IGNORECASE,
        )

        outdir_match = re.search(
            r"\boutdir\s*=\s*['\"]([^'\"]+)['\"]",
            scf_content,
            re.IGNORECASE,
        )

        if prefix_match is None:
            prefix = "pwscf"
        else:
            prefix = prefix_match.group(1)

        if outdir_match is None:
            outdir = "."
        else:
            outdir = outdir_match.group(1)

        proj_input = f"""&PROJWFC
    outdir      = '{outdir}'
    prefix      = '{prefix}'
    filpdos     = 'mh_pdos'
    filproj     = 'mh_proj'
    lsym        = .true.
/
"""

        proj_path = Path(cwd) / "proj.in"

        with open(proj_path, "w") as proj_file:
            proj_file.write(proj_input)

    def _run_projwfc(self, cwd: str = "."):
        """
        Run projwfc.x using the generated proj.in file.
        """

        command = (
            f"{shlex.join(self.projwfc_command)} "
            f"< proj.in"
        )

        subprocess.run(
            command,
            shell=True,
            check=True,
            cwd=cwd,
        )
