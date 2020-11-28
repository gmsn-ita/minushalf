"""
    Creates objects in order to read files potential information
    and correct them.
    Supports the following softwares:
    - VASP
"""
from abc import ABC, abstractmethod
from minushalf.corrections.base_software import Software
from .vasp import Vasp


class CorrectionAbstractFactory(ABC):
    """
    Abstract Factory for create instances for
    each supported software.
    """

    @abstractmethod
    def vasp(self, cut: float, amplitude: float, potcar_path: str) -> Software:
        """
        Concrete method for create instance of
        Vienna Abnition Simulation Package.

         Args:
            cut (float): Value of limit radius necessary for DFT -1/2 techinique
            amplitude (float):  Triming function modulator
            potcar_path (str): Path to POTCAR file

         Returns:
            Vasp class
        """


class CorrectionsFactory(CorrectionAbstractFactory):
    """
    Concrete Factory for create instances
    for each supported software.
    """

    def vasp(self, cut: float, amplitude: float = 1.0, potcar_path: str = "POTCAR") -> Software:
        """
        Concrete method for create instance of
        Vienna Abnition Simulation Package.

         Args:
            cut (float): Value of limit radius necessary for DFT -1/2 techinique
            amplitude (float):  Triming function modulator
            potcar_path (str): Path to POTCAR file

         Returns:
            Vasp class
        """
        return Vasp(cut, amplitude, potcar_path)
