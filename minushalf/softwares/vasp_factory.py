"""
    Factory to generate same modules for different softwares

    Supports the following softwares:
    - VASP
"""
import os
from minushalf.interfaces import SoftwaresAbstractFactory
from minushalf.utils import (
    check_procar_exists,
    check_vasprun_exists,
    check_eigenval_exists,
    check_potcar_exists,
    check_outcar_exists,
)
from minushalf.softwares.vasp import (
    Procar,
    Vasprun,
    Eigenvalues as VaspEigenval,
    Potcar,
    VaspRunner,
    Outcar,
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
                filename (str): Name of the vasprun.xml file.
                base_path (str): Path to the folder where the file is located.

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
                filename (str): Name of the vasprun.xml file.
                base_path (str): Path to the folder where the file is located.

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
                filename (str): Name of the PROCAR file
                base_path (str): Path to the folder where the file is located.

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
                filename (str): Name of the PROCAR file
                base_path (str): Path to the folder where the file is located.

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
                filename (str): Name of the PROCAR file.
                base_path (str): Path to the folder where the file is located.

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
                filename (str): Name of the POTCAR file.
                base_path (str): Path to the folder where the file is located.

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
                filename (str): Name of the EIGENVAL file
                base_path (str): Path to the folder where the file is located.

            Returns:
                eigenvalues (dict): dictionary containing the eigenvalues
                for each kpoint and each band
        """
        if base_path:
            filename = os.path.join(base_path, filename)
        eigenval = VaspEigenval(filename)
        return eigenval.eigenvalues

    def get_runner(self, path: str = "vasp", number_of_cores: int = 1):
        """
        Return the class
        that runs VASP
        """
        return VaspRunner(path, number_of_cores)

    @check_outcar_exists
    def get_nearest_neighbor_distance(self,
                                      ion_index: str,
                                      filename: str = "OUTCAR",
                                      base_path: str = None) -> float:
        """
            Args:
                ion_index (str): The index of the ion given by VASP.
                filename (str): Name of the OUTCAR file.
                base_path (str): Path to the folder where the file is located.

            Returns:
                distance (float): The distance of the nearest neighbor.
        """
        if base_path:
            filename = os.path.join(base_path, filename)

        outcar = Outcar(filename)
        return outcar.nearest_neighbor_distance(ion_index)

    @check_outcar_exists
    def get_number_of_equal_neighbors(self,
                                      atoms_map: dict,
                                      symbol: str,
                                      filename: str = "OUTCAR",
                                      base_path: str = None) -> float:
        """
        Given an map that links atoms symbols with it's index
        this function returns the number of neighbors of the atom with
        equal symbol but different indexes.

            Args:
                atoms_map (dict): Map the atoms index to their symbol.
                symbom (str): The symbol of the target atom.

            Returns:
                number_equal_neighbors (int): Returns the number of neighbors with
                                        same symbol but different indexes.
        """
        if base_path:
            filename = os.path.join(base_path, filename)

        outcar = Outcar(filename)
        return outcar.number_of_equal_neighbors(atoms_map, symbol)
