"""
Leads with input file (INP.ae)
read by atomic program.
"""
import numpy as np
from pathlib import Path
from itertools import dropwhile


class InputFile:
    """
    Parses input file.
    """

    def __init__(self, input_path: str = "./INP") -> None:
        """
        Args:
            input_path (str): Path to INP file
        """
        self.input_path = input_path
        self._parse_file()

    @property
    def input_path(self) -> str:
        """
        Returns:
            Path to INP file.
        """
        return self._input_path

    @input_path.setter
    def input_path(self, input_path: str) -> None:
        """
        Setter for input path, verifies if the file

        Args:
            input_path (str): Path to INP.ae file required by
            atom program.
        """
        path = Path(input_path)

        if not path.exists():
            raise FileNotFoundError()

        self._input_path = path.__str__()

    def electron_occupation(self, occupation_fraction: float) -> None:
        """
        Add the minus one half eléctron correction
        on INP file. Basically, it subtracts half
        electron on the last layer of the eletronic
        structure.
        """
        for value in reversed(self.valence_orbitals):

            if not np.isclose(value["occupation"][0], 0.0, rtol=1e-04, atol=1e-08, equal_nan=False):
                value["occupation"][0] -= occupation_fraction
                break
        else:
            raise Exception("No eletrons in valence orbitals")

    def to_file(self, filename: str = "./INP") -> None:
        """
        Write INP file.
        """
        with open(filename, "w") as file:
            file.write(
                "   {}      {}\n".format(
                    self.calculation_code, self.description)
            )
            file.write(
                " {}".format(self.chemical_simbol)
            )
            if len(self.chemical_simbol) == 4:
                file.write(" {}\n".format(self.exchange_correlation_type))
            else:

                file.write("  {}\n".format(self.exchange_correlation_type))

            file.write("{}".format(self.esoteric_line))

            file.write(
                "    {}    {}\n".format(
                    self.number_core_orbitals, self.number_valence_orbitals
                )
            )
            for orbital in self.valence_orbitals:
                occupation = "      ".join(
                    ["{:.2f}".format(value) for value in orbital["occupation"]]
                )
                file.write(
                    "    {}    {}      {}\n".format(
                        orbital["n"], orbital["l"], occupation
                    )
                )

            for line in self.last_lines:
                file.write("{}".format(line))

    def _parse_file(self):
        """
        Parse INP.ae file.
        """
        with open(self.input_path) as file:
            clean_lines = list(dropwhile(self.is_comment, file))

            if len(clean_lines) < 5:
                raise Exception("Your file is probably incomplete.")

            self.calculation_code = clean_lines[0].split()[0]
            self.description = " ".join(clean_lines[0].split()[1:])
            self.chemical_simbol = clean_lines[1].split()[0]
            self.exchange_correlation_type = clean_lines[1].split()[1]
            self.esoteric_line = clean_lines[2]
            self.number_core_orbitals = int(clean_lines[3].split()[0])
            self.number_valence_orbitals = int(clean_lines[3].split()[1])

            if len(clean_lines) < 4 + self.number_valence_orbitals:
                raise Exception("Your file is probably incomplete.")

            self.valence_orbitals = [
                self._parse_valence_orbitals(clean_lines[i])
                for i in range(4, 4 + self.number_valence_orbitals)
            ]
            self.last_lines = clean_lines[4 + self.number_valence_orbitals:]

    def _parse_valence_orbitals(self, line: str) -> any:
        """Parse valence orbital line in principal
        quantum number, angular momentum quantum number
        and eletronic occupation

        Args:
            line (str): line of imp file that represents
            a valence orbital
        Returns:
            A dictionary with fields n, l and eletronic occupation
        """
        parsed_line = line.split()
        orbital = {
            "n": int(parsed_line[0]),
            "l": int(parsed_line[1]),
            "occupation": [float(value) for value in parsed_line[2:]],
        }

        return orbital

    def is_comment(self, line: str) -> bool:
        """function to check if a line
        starts with some character.
        Here # for comment
        """
        return line.startswith("#")
