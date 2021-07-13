"""
Implementation for vasp runner
"""
import subprocess
from typing import List
from minushalf.interfaces import Runner


class VaspRunner(Runner):
    """
    Output terminal command that
    aims to run vasp
    """
    def __init__(self, command: List[str]):
        """
        Args:
            path_to_vasp (str): Path to vasp binary
            number of cores (int): number of cores needed by the proccess
        """
        self.command = command

    def run(self, cwd: str = "."):
        """
        Create a subproccess to run
        vasp
        """
        subprocess.run(self.command, check=True, cwd=cwd)
