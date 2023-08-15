"""
Band structure informations
"""
from math import inf
from collections import defaultdict
import numpy as np
from minushalf.softwares.band_projection_file import BandProjectionFile
from minushalf.softwares.software_abstract_factory import SoftwaresAbstractFactory



class BandStructure():
    """
    Extact band structure insights from
    VASP classes
    """
    def __init__(self, eigenvalues: dict, fermi_energy: float, atoms_map: dict,
                 num_bands: int, band_projection: BandProjectionFile):
        """
            Args:
                eigenvalues (dict): Egenvalues of dft calculation 
                fermi_energy (float): Value of the fermi level
                atoms_map (dict): List of the atoms of the crystal
                num_bands (int): Number of bands used in the calculation
                band_projection(BandProjectionFile): Class with band projection informations

        """
        self.eigenvalues = eigenvalues
        self.fermi_energy = fermi_energy
        self.atoms_map = atoms_map
        self.num_bands = num_bands
        self.band_projection_file = band_projection

    def is_metal(self, tolerance: float = 1e-4) -> bool:
        """
        Check if the band structure indicates a metal by looking if the fermi
        level crosses a band.
            
            Returns:
                True if a metal, False if not
        """
        eigevalues = [*self.eigenvalues.values()]
        for index in range(self.num_bands):
            band_eigenvalues = np.array([array[index] for array in eigevalues],
                                        dtype=np.float64)
            if np.any(band_eigenvalues -
                      self.fermi_energy < -tolerance) and np.any(
                          band_eigenvalues - self.fermi_energy > tolerance):
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
            raise Exception("Conduction band is not defined for metals")

        kpoint_vbm = None
        band_vbm = None
        max_energy_reached = -inf
        for kpoint, values in self.eigenvalues.items():
            for band_index, energy in enumerate(values):
                if max_energy_reached < energy < self.fermi_energy:
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
            raise Exception("Conduction band is not defined for metals")

        kpoint_cbm = None
        band_cbm = None
        min_energy_reached = inf
        for kpoint, values in self.eigenvalues.items():
            for band_index, energy in enumerate(values):
                if self.fermi_energy <= energy < min_energy_reached:
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
        procar_projection = self.band_projection_file.get_band_projection(
            *vbm_index)

        vbm_projection = {}
        for index, symbol in self.atoms_map.items():
            if vbm_projection.get(symbol, None):
                add_projections = np.add(procar_projection[index],
                                         vbm_projection[symbol])
                vbm_projection[symbol] = list(add_projections)
            else:
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

        procar_projection = self.band_projection_file.get_band_projection(
            *cbm_index)

        cbm_projection = {}
        for index, symbol in self.atoms_map.items():

            if cbm_projection.get(symbol, None):
                add_projections = np.add(procar_projection[index],
                                         cbm_projection[symbol])
                cbm_projection[symbol] = list(add_projections)
            else:
                cbm_projection[symbol] = procar_projection[index]

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
        procar_projection = self.band_projection_file.get_band_projection(
            kpoint, band)

        band_projection = {}

        for index, symbol in self.atoms_map.items():
            if band_projection.get(symbol, None):
                add_projections = np.add(procar_projection[index],
                                         band_projection[symbol])
                band_projection[symbol] = list(add_projections)
            else:
                band_projection[symbol] = procar_projection[index]

        return band_projection

    def band_gap(self) -> dict:
        """
        Find VBM and CBM, then returns band gap
        Returns:
            VBM index and its eigenvalue, CBM index and its eigenvalue and band gap
        """
        vbm = self.vbm_index()
        vbm_eigenval = self.eigenvalues[vbm[0]][vbm[1] - 1]

        cbm = self.cbm_index()
        cbm_eigenval = self.eigenvalues[cbm[0]][cbm[1] - 1]

        gap_report = {
            "vbm":
            "VBM: Kpoint {}, band {} and eigenval {}".format(
                vbm[0], vbm[1], vbm_eigenval),
            "cbm":
            "CBM: Kpoint {}, band {} and eigenval {}".format(
                cbm[0], cbm[1], cbm_eigenval),
            "gap":
            cbm_eigenval - vbm_eigenval
        }
        return gap_report

    @staticmethod
    def create(software_module: SoftwaresAbstractFactory,
               base_path: str = '.'):
        """
        Create band structure class from ab inition results

            Args:
                
                software_module (SoftwaresAbstractFactory): Holds the results of first principles output calculations
                base_path (str): Path to first principles output files
            
            Returns:
                
                band_strucure (BandStructure): Class with band structure informations
        
        """
        eigenvalues = software_module.get_eigenvalues(base_path=base_path)
        fermi_energy = software_module.get_fermi_energy(base_path=base_path)
        atoms_map = software_module.get_atoms_map(base_path=base_path)
        num_bands = software_module.get_number_of_bands(base_path=base_path)
        band_projection_file = software_module.get_band_projection_class(
            base_path=base_path)

        return BandStructure(eigenvalues, fermi_energy, atoms_map, num_bands,
                             band_projection_file)
