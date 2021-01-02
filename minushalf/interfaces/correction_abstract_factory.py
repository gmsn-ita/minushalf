"""
Correction Abstract Factory
"""
from abc import ABC, abstractmethod
from .correction import Correction


class CorrectionAbstractFactory(ABC):
    """
    Abstract factory to generate
    corrections algorithms for
    different softwares
    """
    @abstractmethod
    def valence(self, **kwargs) -> Correction:
        """
        Return a class that does
        valence correction
        """

    @abstractmethod
    def valence_fractionary(self, **kwargs) -> Correction:
        """
        Return a class that does
        fractionary valence correction
        """

    @abstractmethod
    def conduction(self, **kwargs) -> Correction:
        """
        Return a class that does
        conduction correction
        """

    @abstractmethod
    def conduction_fractionary(self, **kwargs) -> Correction:
        """
        Return a class that does
        fractionary conduction correction
        """
