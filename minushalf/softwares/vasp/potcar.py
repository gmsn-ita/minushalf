"""
Reads and analyze POTCAR file
"""
import re
from itertools import chain
import fortranformat as ff
import numpy as np
from minushalf.interfaces.potential_file import PotentialFile


class Potcar(PotentialFile):
    """
    Parse the POTCAR file, a vasp input file. It store
    the fourier coefficients and the restant lines of the file
    """
    def __init__(self, filename: str = "POTCAR") -> None:
        self.filename = filename
        self.psctr_parameters = self._get_psctr_parameters()
        self.k_max_text, self.potential = self._get_potential()
        self.k_max = float(self.k_max_text)
        self.last_lines = self._get_last_lines()
        self.name = "POTCAR"

    def get_potential_fourier_transform(self) -> list:
        """
        Returns:
            potential(list): List of fourier transform
            of the potential
        """
        return self.potential

    def get_name(self) -> str:
        """
        Returns potential file name
        """
        return self.name

    def get_maximum_module_wave_vector(self) -> float:
        """
        Returns:
            k_max (float):maximum modulus of the wave vector in reciprocal space
        """
        return self.k_max

    def to_stringlist(self) -> list:
        """
            Returns:
                potcar_lines (list): List of the POTCAR lines
        """
        fourier_coefficients_lines = [
            "{:<3}{}{:<5}\n".format('', self.k_max_text, '')
        ]
        try:
            grouped_coefficients = self.potential.reshape((-1, 5))
        except ValueError as bad_shape:
            raise ValueError(
                "Fourier coefficents not provided correctly") from bad_shape

        fortran_formater = ff.FortranRecordWriter('(1E26.8)')
        for group in grouped_coefficients:
            formated_numbers = []
            for number in group:
                formated_numbers.append(
                    fortran_formater.write([number]).strip())

            line = "{:<2}{}\n".format('', "  ".join(formated_numbers))
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

    def _get_potential(self) -> tuple:
        """
        Read the maximum fourier coefficient and the list
        of fourier coefficients present in POTCAR

            Returns:
                (k_max,potential) (tuple): A tuple
                containing the maximum modulus of the wave vector in
                reciprocal space string and a list with all fourier
                transforms for the potential.
        """
        k_max_text = ""
        potential = []
        regex_begin_catch = re.compile(r"^\s*local\s+part")
        regex_end_catch = re.compile(r"^\s*[a-z]+")

        with open(self.filename, "r") as potcar:
            ## Jump lines
            for line in potcar:
                if regex_begin_catch.match(line):
                    break

            k_max_text = potcar.readline().strip()

            for line in potcar:
                if regex_end_catch.match(line):
                    potential = list(chain.from_iterable(potential))
                    return (k_max_text, np.array(potential, dtype=np.float64))
                potential.append(line.split())

    def _get_last_lines(self) -> list:
        """
        Read the last lines in POTCAR

            Returns:
                last_lines (list): A list with the lines that are after
                the fourier coefficients. These lines will not be used
                for correction purposes
        """
        last_lines = []
        regex_begin_potential = re.compile(r"^\s*local\s+part")
        regex_end_potential = re.compile(r"^\s*[a-z]+")

        with open(self.filename, "r") as potcar:
            ## Jump potential lines
            for line in potcar:
                if regex_begin_potential.match(line):
                    break

            for line in potcar:
                if regex_end_potential.match(line):
                    last_lines.append(line)
                    break
            ## Catch all lines from here
            for line in potcar:
                last_lines.append(line)

        return last_lines
