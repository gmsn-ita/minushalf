"""
    Factory to generate same modules for different softwares

    Supports the following softwares:
    - VASP
"""
from abc import ABC, abstractmethod
from .vasp import Procar, Vasprun, BandStructure, Eigenvalues


class SoftwaresAbstractFactory(ABC):
    """
    Abstract Factory for create instances for
    each supported software.
    """
    @abstractmethod
    def band_structure(self) -> any:
        """
        Abstract method for create instance of
        band structure class.
        """


class VaspFactory(SoftwaresAbstractFactory):
    """
    Concrete Factory for create instances
    for each supported software.
    """
    def band_structure(self, vasprun_path: str, procar_path: str,
                       eigenval_path: str) -> any:
        """
        Concrete method for create instance band structure function

        Args:
            vasprun_path (str): Path to vasprun.xml
            eigenval_path (str): Path to EIGENVAL
            procar_path (str): Path to PROCAR

        """
        procar = Procar(procar_path)
        vasprun = Vasprun(vasprun_path)
        eigenvalues = Eigenvalues(eigenval_path)

        return BandStructure(procar, vasprun, eigenvalues)
