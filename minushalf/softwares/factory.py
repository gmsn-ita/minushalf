"""
    Factory to generate same modules for different softwares

    Supports the following softwares:
    - VASP
"""
from abc import ABC, abstractmethod
from .vasp import (
    BandStructure as VaspBs,
    Procar,
    Vasprun,
    Eigenvalues as VaspEigenval,
    Potcar,
    AtomicPotential as VaspAtomPotential,
)
from minushalf.atomic import Vtotal, InputFile


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

    @abstractmethod
    def atomic_potential(self) -> any:
        """
        Abstract method for create instance of
        atomic potential class.
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
        eigenvalues = VaspEigenval(eigenval_path)

        return VaspBs(procar, vasprun, eigenvalues)

    def atomic_potential(self, vtotal: Vtotal, vtotal_occupied: Vtotal,
                         input_file: InputFile, potcar_path: str):
        """
        Concrete method for create instance atomic potential function

        Args:
            vtotal (Vtotal): atomic pseudopotential of the atom
            vtotal_occupied (Vtotal): atomic pseudopotential of the occupied atom
            input_file (Vtotal): input file of the occupied atom
            potcar (Potcar): input for VASP
        """
        potcar = Potcar(potcar_path)
        return VaspAtomPotential(vtotal, vtotal_occupied, input_file, potcar)
