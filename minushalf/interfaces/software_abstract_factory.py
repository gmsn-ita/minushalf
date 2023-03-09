"""
Software abstract Factory
"""
from abc import ABC, abstractmethod
from typing import List
from minushalf.interfaces.potential_file import PotentialFile
from minushalf.interfaces.band_projection_file import BandProjectionFile


class SoftwaresAbstractFactory(ABC):
    """
    Abstract Factory for create instances for
    each supported software.
    """
    @abstractmethod
    def get_atoms_map(self, filename: str, base_path: str = None) -> dict:
        """
        Abstract method for returns a map
        of the atomic symbol to its index.
        """

    @abstractmethod
    def get_band_projection_class(
        self,
        filename: str,
        base_path: str = None,
    ) -> BandProjectionFile:
        """
        Abstract method for returns the class that
        handles with the projections of atoms orbitals
        in the bands.
        """

    @abstractmethod
    def get_fermi_energy(self, filename: str, base_path: str = None) -> float:
        """
        Abstract method for returna
        energy of the fermi level.
        """

    @abstractmethod
    def get_number_of_bands(self, filename: str, base_path: str = None) -> int:
        """
        Abstract method for returns the number of bands used in the calculation
        """

    @abstractmethod
    def get_number_of_kpoints(self,
                              filename: str,
                              base_path: str = None) -> int:
        """
        Abstract method for returns the number of kpoints used in the calculation
        """

    @abstractmethod
    def get_potential_class(self,
                            filename: str,
                            base_path: str = None) -> PotentialFile:
        """
        Abstract method for returns the potential class
        """

    @abstractmethod
    def get_eigenvalues(self, filename: str, base_path: str = None) -> dict:
        """
        Abstract method for returns eigenvalues
        for each band and each kpoint
        """

    @abstractmethod
    def get_runner(self, command: List[str]):
        """
        Get a class that run
        the software for ab initio calculations
        """

    @abstractmethod
    def get_nearest_neighbor_distance(self,
                                      ion_index: str,
                                      filename: str,
                                      base_path: str = None) -> float:
        """
        Abstract method for returns the nearest neighbor distance for
        an ion in the solid.
        """

    @abstractmethod
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
