"""
Band structure informations
"""
import re
import numpy as np
from collections import defaultdict
from math import inf


class BandStructure():
    """
    Extact band structure insights from
    VASP classes
    """
    def __init__(self, procar, vasprun, eigenval):
        """
            Args:
                filename (str): name of the EIGENVAL file in VASP
            Members:
        """
        self.procar = procar
        self.vasprun = vasprun
        self.eigenval = eigenval

    def is_metal(self, tolerance: float = 1e-4) -> bool:
        """
        Check if the band structure indicates a metal by looking if the fermi
        level crosses a band.
        Returns:
            True if a metal, False if not
        """
        eigevalues = [*self.eigenval.eigenvalues.values()]
        for index in range(self.procar.num_bands):
            band_eigenvalues = np.array([array[index] for array in eigevalues],
                                        dtype=np.float)
            if np.any(band_eigenvalues -
                      self.vasprun.fermi_energy < -tolerance) and np.any(
                          band_eigenvalues -
                          self.vasprun.fermi_energy > tolerance):
                return True
        return False

    def vbm_index(self) -> tuple:
        """
        Find the kpoint and the band for vbm
            
            Returns:
                vbm_index (tuple): Contains the kpoint
                number and the band number of th vbm
        """
        if self.is_metal():
            raise Error("Conduction band is not defined for metals")

        kpoint_vbm = None
        band_vbm = None
        max_energy_reached = -inf
        for kpoint, values in self.eigenval.eigenvalues.items():
            for band_index, energy in enumerate(values):
                if energy < self.vasprun.fermi_energy and energy > max_energy_reached:
                    max_energy_reached = energy
                    kpoint_vbm = kpoint
                    band_vbm = band_index + 1

        return (kpoint_vbm, band_vbm)

    def cbm_index(self) -> tuple:
        """
        Find the kpoint and the band for cbm 

            Returns:
                vbm_index (tuple): Contains the kpoint
                number and the band number of th vbm
        """
        if self.is_metal():
            raise Error("Conduction band is not defined for metals")

        kpoint_cbm = None
        band_cbm = None
        min_energy_reached = inf
        for kpoint, values in self.eigenval.eigenvalues.items():
            for band_index, energy in enumerate(values):
                if energy >= self.vasprun.fermi_energy and energy < min_energy_reached:
                    min_energy_reached = energy
                    kpoint_cbm = kpoint
                    band_cbm = band_index + 1
        return (kpoint_cbm, band_cbm)

    def vbm_projection(self) -> defaultdict(list):
        """
        Find the projection of each atom for valence
        band maximum.
            
            Returns:
                vbm_projection (defaultdict(list)): Contains the projection
                of each orbital of each atom in the respective band
        """
        vbm_index = self.vbm_index()

        procar_projection = self.procar.get_band_projection(*vbm_index)

        vbm_projection = {}
        for index, symbol in self.vasprun.atoms_map.items():
            vbm_projection[symbol] = procar_projection[index]
        return vbm_projection

    def cbm_projection(self) -> defaultdict(list):
        """
        Find the projection of each atom for valence
        band minimum.
            
            Returns:
                vbm_projection (defaultdict(list)): Contains the projection
                of each orbital of each atom in the respective band
        """
        cbm_index = self.cbm_index()
        procar_projection = self.procar.get_band_projection(*cbm_index)

        cbm_projection = {}
        for index, symbol in self.vasprun.atoms_map.items():
            cbm_projection[symbol] = procar_projection[index]
        return cbm_projection

    def band_projection(self, kpoint: int, band: int) -> defaultdict(list):
        """
        Find the projection of each atom for a specific band.
            Args:
                kpoint (int): Number of kpoints 
                band_number (int): Number of the band
            
            Returns:
                band_projection (defaultdict(list)): Contains the projection
                of each orbital of each atom in the respective band
        """
        procar_projection = self.procar.get_band_projection(kpoint, band)

        band_projection = {}
        for index, symbol in self.vasprun.atoms_map.items():
            band_projection[symbol] = procar_projection[index]
        return band_projection
