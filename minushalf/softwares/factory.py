"""
    Factory to generate same modules for different softwares

    Supports the following softwares:
    - VASP
"""
import os
from abc import ABC, abstractmethod
from minushalf.atomic import Vtotal, InputFile
from .vasp import (
    BandStructure as VaspBs,
    Procar,
    Vasprun,
    Eigenvalues as VaspEigenval,
    Potcar,
    AtomicPotential as VaspAtomPotential,
)


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
    def band_structure(
        self,
        vasprun_path: str = "vasprun.xml",
        procar_path: str = "PROCAR",
        eigenval_path: str = "EIGENVAL",
        base_path: str = None,
    ) -> any:
        """
        Concrete method for create instance band structure function

        Args:
            vasprun_path (str): Path to vasprun.xml
            eigenval_path (str): Path to EIGENVAL
            procar_path (str): Path to PROCAR
            base_path(str): Base path to be apllied in all remaining filenames

        """
        if base_path:
            procar_path = os.path.join(base_path, procar_path)
            vasprun_path = os.path.join(base_path, vasprun_path)
            eigenval_path = os.path.join(base_path, eigenval_path)
        procar = Procar(procar_path)
        vasprun = Vasprun(vasprun_path)
        eigenvalues = VaspEigenval(eigenval_path)

        return VaspBs(procar, vasprun, eigenvalues)

    def atomic_potential(self,
                         vtotal: Vtotal,
                         vtotal_occupied: Vtotal,
                         input_file: InputFile,
                         potcar_path: str = 'POTCAR',
                         base_path: str = None):
        """
        Concrete method for create instance atomic potential function

        Args:
            vtotal (Vtotal): atomic pseudopotential of the atom
            vtotal_occupied (Vtotal): atomic pseudopotential of the occupied atom
            input_file (Vtotal): input file of the occupied atom
            potcar (Potcar): input for VASP
            base_path (str): base_path to potcar
        """
        if base_path:
            potcar_path = os.path.join(base_path, potcar_path)
        potcar = Potcar(potcar_path)
        return VaspAtomPotential(vtotal, vtotal_occupied, input_file, potcar)
