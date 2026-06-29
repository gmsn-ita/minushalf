"""
Reads the pw.x standard output file,
an output of Quantum ESPRESSO software
"""
import re
import numpy as np
from collections import Counter, defaultdict


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
            atoms_map:          dict mapping each atom index (str) to its
                                chemical symbol.
            num_kpoints:        number of k-points used in the simulation
            num_bands:          number of Kohn-Sham states (bands)
            relative_distances: defaultdict(list) mapping each ion index (str)
                                to a list of (neighbor_index, distance_Å) tuples.
        """
        self.filename = filename
        self.atoms_map = self._get_atoms_map()
        self.num_kpoints, self.num_bands = self._get_num_kpoints_and_bands()
        self.relative_distances = self._get_distances()

    # ------------------------------------------------------------------ #
    #  Public neighbor methods                                           #
    # ------------------------------------------------------------------ #

    def nearest_neighbor_distance(self, ion_index: str) -> float:
        """
        Given the ion index, returns the distance of the nearest neighbor
        to this ion in Angstrom.

        Args:
            ion_index (str): 1-based ion index as a string.
        Returns:
            nearest_neighbor_distance (float): distance in Angstrom.
        """
        distances = [dist for _, dist in self.relative_distances[ion_index]]
        return min(distances)

    def number_of_equal_neighbors(self, atoms_map: dict, symbol: str) -> int:
        """
        Given a map that links atom indices to their symbols, returns the
        number of neighbors of the atom with the given symbol that share
        the same symbol but have different indices.

        Args:
            atoms_map (dict): maps atom index (str) to chemical symbol (str).
            symbol    (str):  chemical symbol of the target atom.
        Returns:
            number_equal_neighbors (int)
        """
        ion_index = self._get_ion_index(atoms_map, symbol)
        number_equal_neighbors = 0
        visited_neighbors = defaultdict(bool)

        for index, _ in self.relative_distances[ion_index]:
            index = str(index)
            if self._is_neighbor_equal(index, ion_index,
                                       visited_neighbors, atoms_map, symbol):
                visited_neighbors[index] = True
                number_equal_neighbors += 1

        return number_equal_neighbors

    # ------------------------------------------------------------------ #
    #  Parsers                                                             #
    # ------------------------------------------------------------------ #

    def _get_atoms_map(self) -> dict:
        """
        Extract atom indices and their chemical symbols from the
        'Cartesian axes' block of the pw.x output file.

        Returns:
            atoms_map (dict): { index_str: symbol }
                e.g. { '1': 'Al', '2': 'Al', '3': 'N', '4': 'N' }

        """
        block_lines = self._catch_cartesian_axes_block()

        atom_line_regex = re.compile(
            r'^\s*(\d+)\s+([A-Za-z]+)\s+tau\s*\('
        )

        atoms_map = {}
        for line in block_lines:
            match = atom_line_regex.match(line)
            if match:
                atoms_map[match.group(1)] = match.group(2)

        if not atoms_map:
            raise Exception(
                "PWOutput parser could not find any atoms "
                "in the 'Cartesian axes' block."
            )

        return atoms_map

    def _get_num_kpoints_and_bands(self) -> tuple:
        """
        Extract the number of k-points and number of bands (Kohn-Sham
        states) from the pw.x output file in a single pass.

        Returns:
            (num_kpoints, num_bands) (tuple[int, int])
        """
        kpoints_regex = re.compile(r"^\s*number of k points=\s*([0-9]+)")
        bands_regex   = re.compile(r"^\s*number of Kohn-Sham states\s*=\s*([0-9]+)")

        number_kpoints = None
        number_bands   = None

        with open(self.filename, "r") as fh:
            for line in fh:
                if number_kpoints is None:
                    m = kpoints_regex.match(line)
                    if m:
                        number_kpoints = int(m.group(1))
                if number_bands is None:
                    m = bands_regex.match(line)
                    if m:
                        number_bands = int(m.group(1))
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

    def _get_alat(self) -> float:
        """
        Extract celldm(1) — the lattice parameter in Bohr (alat) — from
        the pw.x output file.

        Returns:
            alat (float): lattice parameter in Bohr
        """
        alat_regex = re.compile(r"^\s*celldm\(1\)=\s*([-+]?[0-9]*\.?[0-9]+)")

        with open(self.filename, "r") as fh:
            for line in fh:
                m = alat_regex.match(line)
                if m:
                    return float(m.group(1))

        raise Exception(
            f"PWOutput parser could not find 'celldm(1)' in {self.filename}"
        )

    def _get_lattice_vectors(self) -> np.ndarray:
        """
        Extract the three crystal (lattice) vectors from the
        'crystal axes' block in units of alat.

        Returns:
            lattice (np.ndarray): shape (3, 3), rows are a1, a2, a3
                                  in units of alat
        """
        header_regex = re.compile(r"^\s*crystal axes:")
        vector_regex = re.compile(
            r"^\s*a\([123]\)\s*=\s*\(\s*"
            r"([-+]?[0-9]*\.?[0-9]+)\s+"
            r"([-+]?[0-9]*\.?[0-9]+)\s+"
            r"([-+]?[0-9]*\.?[0-9]+)\s*\)"
        )

        vectors = []
        inside  = False

        with open(self.filename, "r") as fh:
            for line in fh:
                if not inside:
                    if header_regex.match(line):
                        inside = True
                    continue
                m = vector_regex.match(line)
                if m:
                    vectors.append([float(m.group(1)),
                                    float(m.group(2)),
                                    float(m.group(3))])
                if len(vectors) == 3:
                    break

        if len(vectors) != 3:
            raise Exception(
                "PWOutput parser could not find all three crystal axes "
                f"in {self.filename}"
            )

        return np.array(vectors)

    def _get_atom_positions(self) -> dict:
        """
        Extract Cartesian positions (in alat units) for every atom from
        the 'Cartesian axes' block.

        Returns:
            positions (dict): { index_str: np.ndarray([x, y, z]) }
                              coordinates in units of alat
        """
        block_lines = self._catch_cartesian_axes_block()

        pos_regex = re.compile(
            r'^\s*(\d+)\s+[A-Za-z]+\s+tau\s*\(\s*\d+\s*\)\s*=\s*\(\s*'
            r'([-+]?[0-9]*\.?[0-9]+)\s+'
            r'([-+]?[0-9]*\.?[0-9]+)\s+'
            r'([-+]?[0-9]*\.?[0-9]+)\s*\)'
        )

        positions = {}
        for line in block_lines:
            m = pos_regex.match(line)
            if m:
                positions[m.group(1)] = np.array([
                    float(m.group(2)),
                    float(m.group(3)),
                    float(m.group(4)),
                ])

        if not positions:
            raise Exception(
                "PWOutput parser could not find atom positions "
                "in the 'Cartesian axes' block."
            )

        return positions

    def _get_distances(self) -> defaultdict:
        """
        Compute pairwise distances between all atoms under periodic boundary
        conditions, mirroring Outcar._get_distances.

        Algorithm:
            For each pair (i, j) — including i == j for self-images —
            generate all 27 periodic images of atom j by translating it
            by n1*a1 + n2*a2 + n3*a3 for n1,n2,n3 in {-1, 0, 1}.
            Compute the Euclidean distance from atom i to each image.
            Keep the minimum non-zero distance (zero would be the atom
            with itself in the same cell, i.e. i == j and image (0,0,0)).

        All internal arithmetic is in alat units; distances are converted
        to Angstrom at the end to match Outcar convention.

        Returns:
            relative_distances (defaultdict(list)):
                { ion_index_str: [(neighbor_index_str, distance_Å), ...] }
                sorted by ascending distance for each ion.
        """
        alat      = self._get_alat()
        lattice   = self._get_lattice_vectors()          # (3,3), alat units
        positions = self._get_atom_positions()           # {str: ndarray alat}

        # Pre-build all 27 translation vectors in alat units
        offsets = [
            n1 * lattice[0] + n2 * lattice[1] + n3 * lattice[2]
            for n1 in (-1, 0, 1)
            for n2 in (-1, 0, 1)
            for n3 in (-1, 0, 1)
        ]

        relative_distances = defaultdict(list)

        for i, pos_i in positions.items():
            for j, pos_j in positions.items():
                min_dist = None
                for offset in offsets:
                    diff = pos_i - (pos_j + offset)
                    dist = float(np.linalg.norm(diff))
                    # Skip the trivial self-distance (atom i == atom j,
                    # zero-translation image)
                    if dist < 1e-10:
                        continue
                    if min_dist is None or dist < min_dist:
                        min_dist = dist

                if min_dist is not None:
                    # 1 Bohr = 0.529177 Angstrom - QE outputs in Bohr, but minushalf expects in Angstrom.
                    BOHR_TO_ANGSTROM = 0.529177
                    relative_distances[i].append(
                        (j, min_dist * alat * BOHR_TO_ANGSTROM)
                    )

            # Sort by ascending distance so nearest_neighbor_distance
            # is consistent with Outcar behaviour
            relative_distances[i].sort(key=lambda t: t[1])

        return relative_distances

    # ------------------------------------------------------------------ #
    #  Private helpers                                                   #
    # ------------------------------------------------------------------ #

    def _get_ion_index(self, atoms_map: dict, target_symbol: str) -> str:
        """
        Given the symbol, returns the ion index (str).
        """
        for index, symbol in atoms_map.items():
            if target_symbol == symbol:
                return index

    def _is_neighbor_equal(
        self,
        target_index: str,
        source_index: str,
        visited_neighbors: dict,
        atoms_map: dict,
        symbol: str,
    ) -> bool:
        """
        Checks if elements have the same symbol and different indices
        and have not been visited yet.
        """
        different_index = target_index != source_index
        same_symbol     = atoms_map.get(target_index) == symbol
        not_visited     = not visited_neighbors[target_index]
        return different_index and same_symbol and not_visited

    def _catch_cartesian_axes_block(self) -> list:
        """
        Capture all lines belonging to the 'Cartesian axes' site-position
        block, from the header line up to (but not including) the first
        blank line that follows the atom entries.

        Returns:
            block_lines (list[str]): raw lines of the block
        """
        header_regex    = re.compile(r'^\s*Cartesian axes\s*$')
        atom_line_regex = re.compile(r'tau\s*\(')

        block_lines      = []
        inside_block     = False
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
