"""
Correction interface
"""
from abc import ABC, abstractmethod


class Correction(ABC):
    """
    Interface to correction
    algorithms
    """
    @abstractmethod
    def execute(self):
        """
        Execute the correction
        """
