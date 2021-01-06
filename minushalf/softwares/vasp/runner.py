"""
Implementation for vasp runner
"""
import subprocess
from minushalf.interfaces import Runner


class VaspRunner(Runner):
    """
    Output terminal command that
    aims to run vasp
    """
    def __init__(self,
                 path_to_vasp: str = "vasp",
                 number_of_cores: int = None):
        """
        Args:
            path_to_vasp (str): Path to vasp binary
            number of cores (int): number of cores needed by the proccess
        """
        self.path_to_vasp = path_to_vasp
        if not number_of_cores:
            self.number_of_cores = 1
        else:
            self.number_of_cores = number_of_cores

    def run(self):
        """
        Create a subproccess to run
        vasp
        """
        command = [
            "mpirun", "-np",
            str(self.number_of_cores), self.path_to_vasp
        ]

        subprocess.run(command, check=True)
