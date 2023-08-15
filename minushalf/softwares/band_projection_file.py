"""
Interfaces for files thar
contain band projections
"""
from abc import ABC, abstractmethod


class BandProjectionFile(ABC):
    """
    Class to handle files
    containing informations
    about projections of
    the atoms in bands.
    """
    @abstractmethod
    def get_band_projection(self, kpoint: int, band_number: int) -> dict:
        """
        Abstract method for return the
        contribution of each atom orbital in a
        specific band of an specific kpoint.
        """
