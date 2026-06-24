"""
    Factory to generate same modules for different softwares

    - Quantum Espresso
"""

import os
from typing import List

from minushalf.softwares.software_abstract_factory import SoftwaresAbstractFactory

from minushalf.softwares.qe.pwoutput import PWOutput
from minushalf.softwares.qe.pwscf import PWSCF
from minushalf.softwares.qe.projoutput import ProjOutput
from minushalf.softwares.qe.potential import Potential
from minushalf.softwares.qe.runner import QERunner


class QE(SoftwaresAbstractFactory):
    """
    Concrete Factory for create instances
    for Quantum Espresso
    """

    def get_atoms_map(self,
                      filename: str,
                      base_path: str = None) -> dict:
        """
        Args:
            filename (str): Name of the output file.
            base_path (str): Path to the folder where the file is located.
        Returns:
            atoms_map (dict): Map of atomic symbols to their respective indexes.
        """
        if base_path:
            filename = os.path.join(base_path, filename)
        pwout = PWOutput(filename, self.syst)
        return pwout.atoms_map
    
    def get_fermi_energy(self,
                         filename: str = 'pwscf.xml',
                         base_path: str = None) -> float:
        """
            Args:
                filename (str): Name of the 'prefix'.xml file from QE pw.x scf calculation.
                base_path (str): Path to the folder where the file is located.

            Returns:
                fermi_energy (dict): Energy of the fermi level
        """
        if base_path:
            filename = os.path.join(base_path, filename)
        pwscf = PWSCF(filename, self.syst)
        return pwscf.fermi_energy
    
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
                            filename: str,
                            base_path: str = None) -> int:
        """
            Args:
                filename (str): Name of the output file from pw.x.
                base_path (str): Path to the folder where the file is located.

            Returns:
                number_of_bands(int): Number of bands used in calculation
        """
        if base_path:
            filename = os.path.join(base_path, filename)
        pwout = PWOutput(filename)
        return pwout.num_bands

    def get_number_of_kpoints(self,
                              filename: str,
                              base_path: str = None) -> int:
        """
            Args:
                filename (str): Name of the output file from pw.x.
                base_path (str): Path to the folder where the file is located.

            Returns:
                number_of_kpoints(int): Number of kpoints used in calculation
        """
        if base_path:
            filename = os.path.join(base_path, filename)
        pwout = PWOutput(filename)
        return pwout.num_kpoints

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
                        filename: str,
                        base_path: str = None) -> dict:
        """
            Args:
                filename (str): Name of the 'prefix'.xml file from QE pw.x scf calculation.
                base_path (str): Path to the folder where the file is located.

            Returns:
                eigenvalues (dict): dictionary containing the eigenvalues
                for each kpoint and each band
        """
        if base_path:
            filename = os.path.join(base_path, filename)
        pwscf = PWSCF(filename)
        return pwscf.eigenvalues

    def get_runner(self, command: List[str]):
        """
        Return the class
        that runs QE pw.x
        > Missing implementation to run ld1.x and virtual_v2.x
        """
        return QERunner(command)
