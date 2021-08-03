"""
Get band projections
"""
import collections
import loguru
import pandas as pd
from minushalf.softwares import Vasp
from minushalf.data import Softwares
from minushalf.handlers import BaseHandler
from minushalf.interfaces import SoftwaresAbstractFactory, BandProjectionFile, MinushalfYaml
from minushalf.utils import (BandProjectionFile, band_structure, projection_to_df)


class GetBandProjections(BaseHandler):
    """
    Uses the software module to extract the character of the bands (VBM and CBM)
    """
    def _get_band_structure(
            self, software_module: SoftwaresAbstractFactory) -> BandProjectionFile:
        """
        Return band structure class
        """
        eigenvalues = software_module.get_eigenvalues()
        fermi_energy = software_module.get_fermi_energy()
        atoms_map = software_module.get_atoms_map()
        num_bands = software_module.get_number_of_bands()
        band_projection_file = software_module.get_band_projection_class()

        return BandProjectionFile(eigenvalues, fermi_energy, atoms_map, num_bands,
                             band_projection_file)

    def _get_vbm_projection(self,
                            band_structure: BandProjectionFile) -> pd.DataFrame:
        """
        Returns vbm projection
        """
        vbm_projection = band_structure.vbm_projection()
        normalized_df = projection_to_df(vbm_projection)
        return normalized_df

    def _get_cbm_projection(self,
                            band_structure: BandProjectionFile) -> pd.DataFrame:
        """
        Returns cbm projection
        """
        cbm_projection = band_structure.cbm_projection()
        normalized_df = projection_to_df(cbm_projection)
        return normalized_df

    def _get_band_projection(self, band_structure: BandProjectionFile, kpoint: int,
                             band: int) -> pd.DataFrame:
        """
        Returns band projection
        """
        band_projection = band_structure.band_projection(kpoint, band)
        normalized_df = projection_to_df(band_projection)
        return normalized_df

    def _select_vbm(self, minushalf_yaml: MinushalfYaml,
                    band_structure: BandProjectionFile) -> pd.DataFrame:
        """
        Select and returns vbm character
        """
        overwrite_vbm = minushalf_yaml.get_overwrite_vbm()

        if overwrite_vbm:
            return self._get_band_projection(band_structure, *overwrite_vbm)

        return self._get_vbm_projection()

    def _select_cbm(self, minushalf_yaml: MinushalfYaml,
                    band_structure: BandProjectionFile) -> pd.DataFrame:
        """
        Select and returns cbm character
        """
        overwrite_cbm = minushalf_yaml.get_overwrite_cbm()

        if overwrite_cbm:
            return self._get_band_projection(band_structure, *overwrite_cbm)

        return self._get_cbm_projection()

    def _is_valence_correction(self, correction_code: str) -> bool:
        """
        Verify if the correction  is a valence correction
        """
        return "v" in correction_code

    def _is_conduction_correction(self, correction_code: str) -> bool:
        """
        Verify if the correction  is a conduction correction
        """
        return "c" in correction_code

    def _get_projections(
            self, minushalf_yaml: MinushalfYaml,
            band_structure: BandProjectionFile) -> collections.defaultdict:
        """
        Returns an dictionary with the projections necessary to the correction
        """
        correction_code = minushalf_yaml.get_correction_code()

        projections = collections.defaultdict(dict)

        ## Select confuction and valence index
        if self._is_valence_correction(correction_code):
            projections["valence"] = self._select_vbm(minushalf_yaml,
                                                      band_structure)
        if self._is_conduction_correction(correction_code):
            projections["conduction"] = self._select_cbm(
                minushalf_yaml, band_structure)
        return projections

    def action(self, request: any) -> any:
        """
        Uses the software module to extract band projections
        """
        loguru.logger.info("Extracting band projections")
        software_module, minushalf_yaml = request["software_module"], request[
            "minushalf_yaml"]
        band_structure = self._get_band_structure(software_module)

        ## Add projections to object
        request["projections"] = self._get_projections(minushalf_yaml,
                                                       band_structure)

        return request
