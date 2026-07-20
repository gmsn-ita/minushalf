"""
Abstract factory interface for creating software-specific objects used in the
minushalf workflow, such as runners, parsers, and pseudopotential handlers.

Supported software packages:
    - Quantum ESPRESSO
"""

import os
from typing import List

from minushalf.softwares.software_abstract_factory import SoftwaresAbstractFactory

from minushalf.softwares.qe.pwscf import PWSCF
from minushalf.softwares.qe.projoutput import ProjOutput
from minushalf.softwares.qe.potential import Potential
from minushalf.softwares.qe.runner import QERunner


class QE(SoftwaresAbstractFactory):
    """
    Concrete implementation of `SoftwaresAbstractFactory` for Quantum ESPRESSO.

    Provides factory methods to create the objects required to run and parse
    Quantum ESPRESSO calculations within the minushalf workflow, including
    runners, parsers, and pseudopotential handlers.
    """
    def __init__(self):
        self._pwscf_cache = {}

    def _load_pwscf(self, filename='pwscf.xml', base_path=None):
        if base_path:
            filename = os.path.join(base_path, filename)
        if filename not in self._pwscf_cache:
            self._pwscf_cache[filename] = PWSCF(filename)
        return self._pwscf_cache[filename]
    

    def get_atoms_map(self,
                      filename: str = 'pwscf.xml',
                      base_path: str = None) -> dict:
        """
        Args:
            filename (str): Name of the 'prefix'.xml file from QE pw.x scf calculation.
            base_path (str): Path to the folder where the file is located.
        Returns:
            atoms_map (dict): Map of atomic symbols to their respective indexes.
        """

        return self._load_pwscf(filename, base_path).atoms_map
    
    def get_fermi_energy(self,
                         filename: str  = 'pwscf.xml',
                         base_path: str = None) -> float:
        """
            Args:
                filename (str): Name of the 'prefix'.xml file from QE pw.x scf calculation.
                base_path (str): Path to the folder where the file is located.

            Returns:
                fermi_energy (dict): Energy of the fermi level
        """
        return self._load_pwscf(filename, base_path).fermi_energy
    
    def get_band_projection_class(
        self,
        filename: str,
        base_path: str = None,
    ) -> ProjOutput:
        """
            Args:
                filename (str): Name of the output file from projwfc.x.
                base_path (str): Path to the folder where the file is located.

            Returns:
                procar (Procar): Contains the class that handles files
                that contains informations about band projections
        """
        if base_path:
            filename = os.path.join(base_path, filename)
        return ProjOutput(filename)
    
    def get_number_of_bands(self,
                            filename: str = 'pwscf.xml',
                            base_path: str = None) -> int:
        """
            Args:
                filename (str): Name of the 'prefix'.xml file from QE pw.x scf calculation.
                base_path (str): Path to the folder where the file is located.

            Returns:
                number_of_bands(int): Number of bands used in calculation
        """
        return self._load_pwscf(filename, base_path).number_of_bands

    def get_number_of_kpoints(self,
                              filename: str = 'pwscf.xml',
                              base_path: str = None) -> int:
        """
            Args:
                filename (str): Name of the 'prefix'.xml file from QE pw.x scf calculation.
                base_path (str): Path to the folder where the file is located.

            Returns:
                number_of_kpoints(int): Number of kpoints used in calculation
        """
        return self._load_pwscf(filename, base_path).number_of_kpoints

    def get_potential_class(
        self,
        filename: str,
        base_path: str = None,
    ) -> Potential:
        """
            Args:
                filename (str): Name of the potential <element>.upf file.
                base_path (str): Path to the folder where the file is located.

            Returns:
                Potcar: class to the potential file
        """
        if base_path:
            filename = os.path.join(base_path, filename)
        return Potential(filename)

    def get_eigenvalues(self,
                        filename: str = 'pwscf.xml',
                        base_path: str = None) -> dict:
        """
            Args:
                filename (str): Name of the 'prefix'.xml file from QE pw.x scf calculation.
                base_path (str): Path to the folder where the file is located.

            Returns:
                eigenvalues (dict): dictionary containing the eigenvalues
                for each kpoint and each band
        """
        return self._load_pwscf(filename, base_path).eigenvalues

    def get_runner(self, command: List[str]):
        """
        Return the class
        that runs QE pw.x
        > Missing implementation to run ld1.x and virtual_v2.x
        """
        return QERunner(command)
    
    
    def get_nearest_neighbor_distance(self,
                                      ion_index: str,
                                      filename: str = 'pwscf.xml',
                                      base_path: str = None) -> float:
        """
            Args:
                ion_index (str): The index of the ion from atoms map.
                filename (str): Name of the 'prefix'.xml file from QE pw.x scf calculation.
                base_path (str): Path to the folder where the file is located.

            Returns:
                distance (float): The distance of the nearest neighbor.
        """
        return self._load_pwscf(filename, base_path).nearest_neighbor_distance(ion_index)

    def get_number_of_equal_neighbors(self,
                                      atoms_map: dict,
                                      symbol: str,
                                      filename: str = 'pwscf.xml',
                                      base_path: str = None) -> float:
        """
        Given an map that links atoms symbols with it's index
        this function returns the number of neighbors of the atom with
        equal symbol but different indexes.

            Args:
                atoms_map (dict): Map the atoms index to their symbol.
                symbom (str): The symbol of the target atom.
                filename (str): Name of the 'prefix'.xml file from QE pw.x scf calculation.
                base_path (str): Path to the folder where the file is located.

            Returns:
                number_equal_neighbors (int): Returns the number of neighbors with
                                        same symbol but different indexes.
        """
        return self._load_pwscf(filename, base_path).number_of_equal_neighbors(atoms_map, symbol)

