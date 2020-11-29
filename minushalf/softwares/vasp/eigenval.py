"""
Reads Eigenval file, an output of
VASP software
"""
import re
from collections import defaultdict


class Eigenvalues():
    """
    Reads eigenvalues and store
    useful informations
    """
    def __init__(self, filename: str):
        """
            Args:
                filename (str): name of the EIGENVAL file in VASP
            Members:
                eigenvalues (defaultdict(list)): An dictionary where
                the keys are kpoint index and the values are lists
                with the eigenvalues for each band
        """
        self.filename = filename
        self.eigenvalues = self._get_eigenvalues()

    def _get_eigenvalues(self) -> defaultdict(list):
        """
            Returns:
                eigenvalues (defaultdict(list)): An dictionary where
                the keys are kpoint index and the values are lists
                with the eigenvalues for each band
        """
        eigenvalues = defaultdict(list)
        kpoint = 0
        band = None
        eigenvalue_regex = re.compile(r"\s*([0-9]+)\s+([-+]?[0-9]*\.[0-9]+)")
        with open(self.filename) as eigenval:
            for line in eigenval:
                if eigenvalue_regex.match(line):
                    band = eigenvalue_regex.match(line).group(1)
                    eigenvalue = eigenvalue_regex.match(line).group(2)

                    if band == 1:
                        kpoint += 1

                    eigenvalues[kpoint].append(eigenvalue)
        return eigenvalues
