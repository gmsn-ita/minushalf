"""
Reads and analyze POTCAR file
"""
import re
from itertools import chain
import numpy as np


class Potcar():
    """
    Parse the POTCAR file, a vasp input file. It store
    the fourier coefficients and the restant lines of the file
    """
    def __init__(self, filename: str = "POTCAR") -> None:
        self.filename = filename
        self.psctr_parameters = self._get_psctr_parameters()
        self.k_max, self.fourier_coefficients = self._get_fourier_coefficients(
        )
        self.last_lines = self._get_last_lines()

    def to_stringlist(self) -> list:
        """
            Returns:
                potcar_lines (list): List of the POTCAR lines
        """
        fourier_coefficients_lines = ["   {:.18f}     ".format(self.k_max)]
        try:
            grouped_coefficients = self.fourier_coefficients.reshape((-1, 5))
        except ValueError as bad_shape:
            raise ValueError(
                "Fourier coefficents not provided correctly") from bad_shape

        for group in grouped_coefficients:
            line = "  {}  {}  {}  {}  {}".format(
                format(group[0], "15.8E"),
                format(group[1], "15.8E"),
                format(group[2], "15.8E"),
                format(group[3], "15.8E"),
                format(group[4], "15.8E"),
            )
            fourier_coefficients_lines.append(line)

        lines = [
            self.psctr_parameters, fourier_coefficients_lines, self.last_lines
        ]
        return list(chain.from_iterable(lines))

    def to_file(self, filename: str) -> None:
        """
        Write POTCAR file
            Args:
                filename (str): Name of the file
        """
        lines = self.to_stringlist()
        with open(filename, "w") as new_potcar:
            new_potcar.writelines(lines)

    def _get_psctr_parameters(self) -> list:
        """
        Read PSCTR parameters in POTCAR

            Returns:
                psctr_parameters (list): A list with the lines that contain
                the psctr_patameters.
        """
        psctr_parameters = []
        regex_end_psctr_catch = re.compile(r"^\s*local\s+part")
        with open(self.filename, "r") as potcar:
            for line in potcar:
                psctr_parameters.append(line)
                if regex_end_psctr_catch.match(line):
                    return psctr_parameters

    def _get_fourier_coefficients(self) -> tuple:
        """
        Read the maximum fourier coefficient and the list
        of fourier coefficients present in POTCAR

            Returns:
                (k_max,fourier coefficients) (tuple): A tuple
                containing the maximum fourier coefficient and
                a list with all cofficients.
        """
        k_max = None
        fourier_coefficients = []
        regex_begin_catch = re.compile(r"^\s*local\s+part")
        regex_end_catch = re.compile(
            r"^\s*gradient\s+corrections\s+used\s+for\s+XC")

        with open(self.filename, "r") as potcar:
            ## Jump lines
            for line in potcar:
                if regex_begin_catch.match(line):
                    break

            k_max = float(potcar.readline())

            for line in potcar:
                if regex_end_catch.match(line):
                    fourier_coefficients = list(
                        chain.from_iterable(fourier_coefficients))
                    return (k_max,
                            np.array(fourier_coefficients, dtype=np.float))
                fourier_coefficients.append(line.split())

    def _get_last_lines(self) -> list:
        """
        Read the last lines in POTCAR

            Returns:
                last_lines (list): A list with the lines that are after
                the fourier coefficients. These lines will not be used
                for correction purposes
        """
        last_lines = []
        regex_begin_catch = re.compile(
            r"^\s*gradient\s+corrections\s+used\s+for\s+XC")
        with open(self.filename, "r") as potcar:
            for line in potcar:
                if regex_begin_catch.match(line):
                    last_lines.append(line)
                    break
            for line in potcar:
                last_lines.append(line)
        return last_lines
