"""
Runner abstract class
"""
from abc import ABC, abstractmethod


class Runner(ABC):
    """
    Run softwares that
    realizes ab initio calculations
    """
    @abstractmethod
    def run(self, cwd: str):
        """
        Create command and
        run the subprocess for it
        """
