"""
Interface for minushalf.yaml tag classes
"""
from abc import ABC, abstractmethod


class MinushalfYamlTags(ABC):
    """
    Interface
    """
    @abstractmethod
    def to_list(self):
        """
        Get list of variables
        """

    @abstractmethod
    def to_dict(self):
        """
        Get dictionary of variables
        """
