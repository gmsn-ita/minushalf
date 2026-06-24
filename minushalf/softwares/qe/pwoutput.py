"""
Reads the pw.x standard output file,
an output of Quantum ESPRESSO software
"""
import re
from collections import Counter


class PWOutput():
    """
    Reads a pw.x standard output file (e.g. pwscf.out) and stores
    information parsed from it.
    """

    def __init__(self, filename: str):
        """
        Args:
            filename (str): path to the pw.x standard output file
        Members:
            atoms_map:   dict mapping each atomic species to its count
                         e.g. { 'Al': 2, 'N': 2 }
            num_kpoints: number of k-points used in the simulation
            num_bands:   number of Kohn-Sham states (bands) used in the simulation
        """
        self.filename = filename
        self.atoms_map = self._get_atoms_map()
        self.num_kpoints, self.num_bands = self._get_num_kpoints_and_bands()

    # ------------------------------------------------------------------ #
    #  Parsers                                                             #
    # ------------------------------------------------------------------ #

    def _get_atoms_map(self) -> dict:
        """
        Extract atomic species and their counts from the
        'Cartesian axes' block of the pw.x output file.

        Returns:
            atoms_map (dict): { symbol (str): count (int) }
                e.g. { 'Al': 2, 'N': 2 }
        """
        block_lines = self._catch_cartesian_axes_block()

        atom_line_regex = re.compile(
            r'^\s*\d+\s+([A-Za-z]+)\s+tau\s*\('
        )

        species = []
        for line in block_lines:
            match = atom_line_regex.match(line)
            if match:
                species.append(match.group(1))

        if not species:
            raise Exception(
                "PWOutput parser could not find any atoms "
                "in the 'Cartesian axes' block."
            )

        return dict(Counter(species))

    def _get_num_kpoints_and_bands(self) -> tuple:
        """
        Extract the number of k-points and the number of bands (Kohn-Sham
        states) from the pw.x output file in a single pass.

        Returns:
            (num_kpoints, num_bands) (tuple[int, int])
        """
        kpoints_regex = re.compile(
            r"^\s*number of k points=\s*([0-9]+)"
        )
        bands_regex = re.compile(
            r"^\s*number of Kohn-Sham states\s*=\s*([0-9]+)"
        )

        number_kpoints = None
        number_bands = None

        with open(self.filename, "r") as pwscf:
            for line in pwscf:
                if number_kpoints is None:
                    kpoints_match = kpoints_regex.match(line)
                    if kpoints_match:
                        number_kpoints = int(kpoints_match.group(1))

                if number_bands is None:
                    bands_match = bands_regex.match(line)
                    if bands_match:
                        number_bands = int(bands_match.group(1))

                if number_kpoints is not None and number_bands is not None:
                    break

        if number_kpoints is None:
            raise Exception(
                "PWOutput parser could not find 'number of k points' "
                f"in {self.filename}"
            )
        if number_bands is None:
            raise Exception(
                "PWOutput parser could not find 'number of Kohn-Sham states' "
                f"in {self.filename}"
            )

        return (number_kpoints, number_bands)

    # ------------------------------------------------------------------ #
    #  Private helpers                                                     #
    # ------------------------------------------------------------------ #

    def _catch_cartesian_axes_block(self) -> list:
        """
        Capture all lines belonging to the 'Cartesian axes' site-position
        block, from the header line up to (but not including) the first
        blank line that follows the atom entries.

        Returns:
            block_lines (list[str]): raw lines of the block
        """
        header_regex = re.compile(r'^\s*Cartesian axes\s*$')
        atom_line_regex = re.compile(r'tau\s*\(')

        block_lines = []
        inside_block = False
        found_first_atom = False

        with open(self.filename, "r") as fh:
            for line in fh:
                if not inside_block:
                    if header_regex.match(line):
                        inside_block = True
                        block_lines.append(line)
                    continue

                block_lines.append(line)

                if atom_line_regex.search(line):
                    found_first_atom = True
                elif found_first_atom and line.strip() == "":
                    break

        if not block_lines:
            raise Exception(
                "PWOutput parser could not locate the "
                "'Cartesian axes' block in the output file."
            )

        return block_lines
    