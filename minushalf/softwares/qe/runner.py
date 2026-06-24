"""
Implementation for quantum espresso runner
"""
import subprocess
from typing import List
from minushalf.softwares.runner import Runner


class QERunner(Runner):
    """
    Output terminal command that
    aims to runs quantum espresso
    """
    def __init__(self, command: List[str]):
        """
        Args:
            command: the command line from the .yaml file
        """
        self.command = command

    def run(self, cwd: str = "."):
        """
        Create a subproccess to run
        quantum espresso
        """
        subprocess.run(self.command, check=True, cwd=cwd)
