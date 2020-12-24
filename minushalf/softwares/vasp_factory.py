"""
    Factory to generate same modules for different softwares

    Supports the following softwares:
    - VASP
"""
import os
from minushalf.interfaces import SoftwaresAbstractFactory
from minushalf.utils import (check_procar_exists, check_vasprun_exists,
                             check_eigenval_exists, check_potcar_exists)
from .vasp import (
    Procar,
    Vasprun,
    Eigenvalues as VaspEigenval,
    Potcar,
)


class VaspFactory(SoftwaresAbstractFactory):
    """
    Concrete Factory for create instances
    for each supported software.
    """
    @check_vasprun_exists
    def get_atoms_map(self,
                      filename: str = "vasprun.xml",
                      base_path: str = None) -> dict:
        """
            Args:
                filename (str): Path to vasprun.xml
                base_path (str): In case you do not need to specify
                the file name, go to the directory where it is located
            Returns:
                atoms_map (dict): Map of atomic symbols for their respective indexes
        """
        if base_path:
            filename = os.path.join(base_path, filename)
        vasprun = Vasprun(filename)
        return vasprun.atoms_map

    @check_vasprun_exists
    def get_fermi_energy(self,
                         filename: str = "vasprun.xml",
                         base_path: str = None) -> float:
        """
        Args:
                filename (str): Path to vasprun.xml
                base_path (str): In case you do not need to specify
                the file name, go to the directory where it is located
            Returns:
                fermi_energy (dict): Energy of the fermi level
        """
        if base_path:
            filename = os.path.join(base_path, filename)
        vasprun = Vasprun(filename)
        return vasprun.fermi_energy

    @check_procar_exists
    def get_band_projection_class(
        self,
        filename: str = "PROCAR",
        base_path: str = None,
    ) -> Procar:
        """
            Args:
                filename (str): Path to PROCAR
                base_path (str): In case you do not need to specify
                the file name, go to the directory where it is located
            Returns:
                procar (Procar): Contains the class that handles files
                that contains informations about band projections
        """
        if base_path:
            filename = os.path.join(base_path, filename)
        return Procar(filename)

    @check_procar_exists
    def get_number_of_bands(self,
                            filename: str = "PROCAR",
                            base_path: str = None) -> int:
        """
            Args:
                filename (str): Path to PROCAR
                base_path (str): In case you do not need to specify
                the file name, go to the directory where it is located
            Returns:
                number_of_bands(int): Number of bands used in calculation
        """
        if base_path:
            filename = os.path.join(base_path, filename)
        procar = Procar(filename)
        return procar.num_bands

    @check_procar_exists
    def get_number_of_kpoints(self,
                              filename: str = "PROCAR",
                              base_path: str = None) -> int:
        """
            Args:
                filename (str): Path to PROCAR
                base_path (str): In case you do not need to specify
                the file name, go to the directory where it is located
            Returns:
                number_of_kpoints(int): Number of kpoints used in calculation
        """
        if base_path:
            filename = os.path.join(base_path, filename)
        procar = Procar(filename)
        return procar.num_kpoints

    @check_potcar_exists
    def get_potential_class(
        self,
        filename: str = "POTCAR",
        base_path: str = None,
    ) -> Potcar:
        """
            Args:
                filename (str): Path to POTCAR
                base_path (str): In case you do not need to specify
                the file name, go to the directory where it is located
            Returns:
                Potcar: class to the potential file
        """
        if base_path:
            filename = os.path.join(base_path, filename)
        return Potcar(filename)

    @check_eigenval_exists
    def get_eigenvalues(self,
                        filename: str = "EIGENVAL",
                        base_path: str = None) -> dict:
        """
            Args:
                filename (str): Path to EIGENVAL
                base_path (str): In case you do not need to specify
                the file name, go to the directory where it is located
            Returns:
                eigenvalues (dict): dictionary containing the eigenvalues
                for each kpoint and each band
        """
        if base_path:
            filename = os.path.join(base_path, filename)
        eigenval = VaspEigenval(filename)
        return eigenval.eigenvalues
