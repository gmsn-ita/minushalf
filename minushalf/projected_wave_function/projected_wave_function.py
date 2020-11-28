"""
    Creates objects in order to read files with static information
    from files containing wave projections.

    Supports the following softwares:
    - VASP
"""
from abc import ABC, abstractmethod
from .vasp import Vasp


class ProjectedWaveFunctionAbstractFactory(ABC):
    """
    Abstract Factory for create instances for
    each supported software.
    """

    @abstractmethod
    def vasp(self, vasprun_path: str) -> Vasp:
        """
        Abstract method for create instance of
        Vienna Abnition Simulation Package.

        Args:
            vasprun_path (str): Path to vasprun.xml generated
            by VASP.
        """


class ProjectedWaveFunctionFactory(ProjectedWaveFunctionAbstractFactory):
    """
    Concrete Factory for create instances
    for each supported software.
    """

    def vasp(self, vasprun_path: str = "./vasprun.xml") -> Vasp:
        """
        Concrete method for create instance of
        Vienna Abnition Simulation Package.

        Args:
            vasprun_path (str): Path to vasprun.xml generated
            by VASP.
        """
        return Vasp(vasprun_path)
