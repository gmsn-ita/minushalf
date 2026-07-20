"""
Implementation of the runner interface for Quantum ESPRESSO.

Provides the concrete runner class responsible for executing pw.x
"""
import subprocess
from typing import List
from minushalf.softwares.runner import Runner


class QERunner(Runner):
    """
    Runner for Quantum ESPRESSO calculations.
    
    Constructs and executes the command-line call to the Quantum ESPRESSO
    binary (e.g. `pw.x`).
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
