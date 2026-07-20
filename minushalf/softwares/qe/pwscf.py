"""
Read the pw.x XML data file ($prefix.xml),
an output of Quantum ESPRESSO software
"""
import numpy as np
import xml.etree.ElementTree as ET
from collections import defaultdict
from typing import Any, Callable

_BOHR_TO_ANGSTROM = 0.529177


class PWSCF():
    """
    Read a pw.x XML output file (e.g. pwscf.xml) and stores
    information parsed from it.
    """

    def __init__(self, filename: str):
        """
        Args:
            filename (str): path to the pw.x XML output file ($prefix.xml)
        Members:
            filename:          stored path
            atoms_map:         dict mapping each atom index (str) to its symbol
                               e.g. { '1': 'Al', '2': 'Al', '3': 'N', '4': 'N' }
            fermi_energy:      Fermi energy in Hartree atomic units
            number_of_kpoints: number of k-points used in the calculation
            number_of_bands:   number of bands used in the calculation
            eigenvalues:       defaultdict(list) where keys are 1-based kpoint
                               indices and values are lists of eigenvalues in
                               Hartree atomic units for each band, in order
            relative_distances: defaultdict(list) mapping each ion index (str)
                               to a list of (neighbor_index, distance_Å) tuples,
                               mirroring Outcar.relative_distances
        """
        self.filename = filename
        self._root = ET.parse(filename).getroot()
        self.atoms_map = self._get_atoms_map()
        self.fermi_energy = self._get_fermi_energy()
        self.number_of_kpoints = self._get_number_of_kpoints()
        self.number_of_bands = self._get_number_of_bands()
        self.eigenvalues = self._get_eigenvalues()
        self.relative_distances = self._get_distances()

    # ------------------------------------------------------------------ #
    #  Parsers                                                             #
    # ------------------------------------------------------------------ #

    def _get_atoms_map(self) -> dict:
        """
        Extract atom indices and their chemical symbols from the
        <atomic_positions> block inside <atomic_structure> of the QE XML file.

        The relevant snippet looks like:
            <atomic_positions>
                <atom name="Al" index="1"> ... </atom>
                <atom name="Al" index="2"> ... </atom>
                <atom name="N"  index="3"> ... </atom>
                <atom name="N"  index="4"> ... </atom>
            </atomic_positions>

        Returns:
            atoms_map (dict): { index_str: symbol }
                e.g. { '1': 'Al', '2': 'Al', '3': 'N', '4': 'N' }
        """
        atoms_map = {}
        for atom in self._root.iter("atom"):
            atoms_map[atom.get("index")] = atom.get("name")

        if not atoms_map:
            raise Exception(
                "PWSCF parser could not find any <atom> elements "
                f"in {self.filename}"
            )
        return atoms_map

    def _get_fermi_energy(self) -> float:
        """
        Extract the Fermi energy in Hartree atomic units from
        <band_structure>/<fermi_energy>.
        """
        return self._get_xml_value(".//fermi_energy", float, "<fermi_energy>")

    def _get_number_of_kpoints(self) -> int:
        """
        Extract the number of k-points from <band_structure>/<nks>.
        """
        return self._get_xml_value(".//band_structure/nks", int,
                                   "<band_structure>/<nks>")

    def _get_number_of_bands(self) -> int:
        """
        Extract the number of bands from <band_structure>/<nbnd>.
        """
        return self._get_xml_value(".//band_structure/nbnd", int,
                                   "<band_structure>/<nbnd>")

    def _get_xml_value(self, xpath: str, cast: Callable[[str], Any],
                       label: str) -> Any:
        """
        Look up a value in the parsed XML tree and cast it to the desired type.

        Args:
            xpath (str): XPath expression locating the element to read.
            cast  (Callable[[str], Any]): callable to convert the element's
                text content to the desired type (e.g. float, int).
            label (str): human-readable label used in error messages.

        Returns:
            The element's text content converted using cast.
        """
        element = self._root.find(xpath)
        if element is None:
            raise Exception(
                f"PWSCF parser could not find {label} in {self.filename}"
            )
        try:
            return cast(element.text)
        except ValueError as invalid_conversion:
            raise Exception(
                f"PWSCF parser could not parse the {label} value"
            ) from invalid_conversion

    def _get_eigenvalues(self) -> defaultdict:
        """
        Extract eigenvalues for every k-point and band from the QE XML file.

        The XML contains one <ks_energies> block per k-point. Inside each
        block, the <eigenvalues> tag holds all band values as space-separated
        floats, in band order.

        Mirrors Eigenvalues._get_eigenvalues from the VASP parser:
          - keys   : 1-based k-point index (int)
          - values : list of floats, one per band in ascending band order,
                     in Hartree atomic units

        Returns:
            eigenvalues (defaultdict(list)):
                { kpoint_index: [e_band1, e_band2, ...], ... }
        """
        eigenvalues = defaultdict(list)
        for kpoint, ks_block in enumerate(self._root.iter("ks_energies"),
                                          start=1):
            eig_element = ks_block.find("eigenvalues")
            if eig_element is not None and eig_element.text:
                eigenvalues[kpoint].extend(
                    float(v) for v in eig_element.text.split()
                )

        if not eigenvalues:
            raise Exception(
                "PWSCF parser could not find any <ks_energies> blocks "
                f"in {self.filename}"
            )
        return eigenvalues

    def _get_distances(self) -> defaultdict:
        """
        Compute pairwise distances between all atoms under periodic boundary
        conditions, using lattice vectors and positions from <atomic_structure>.

        Positions and lattice vectors are in Bohr (Cartesian).
        All 27 periodic images are checked for each pair; the minimum
        non-zero distance is kept and converted to Angstrom.

        Returns:
            relative_distances (defaultdict(list)):
                { ion_index_str: [(neighbor_index_str, distance_Å), ...] }
                sorted by ascending distance for each ion.
        """
        atomic_structure = self._root.find(".//atomic_structure")
        if atomic_structure is None:
            raise Exception(
                f"PWSCF parser could not find <atomic_structure> in {self.filename}"
            )

        # Parse atom positions (Bohr, Cartesian)
        positions = {}
        for atom in atomic_structure.iter("atom"):
            coords = [float(v) for v in atom.text.split()]
            positions[atom.get("index")] = np.array(coords)

        # Parse lattice vectors (Bohr, Cartesian)
        cell = atomic_structure.find("cell")
        lattice = np.array([
            [float(v) for v in cell.find(tag).text.split()]
            for tag in ("a1", "a2", "a3")
        ])

        # Pre-build all 27 translation vectors
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
                    dist = float(np.linalg.norm(pos_i - (pos_j + offset)))
                    if dist < 1e-10:
                        continue
                    if min_dist is None or dist < min_dist:
                        min_dist = dist
                if min_dist is not None:
                    relative_distances[i].append(
                        (j, min_dist * _BOHR_TO_ANGSTROM)
                    )
            relative_distances[i].sort(key=lambda t: t[1])

        return relative_distances

    # ------------------------------------------------------------------ #
    #  Public neighbor methods  (mirrors Outcar interface)                 #
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
    #  Private helpers                                                     #
    # ------------------------------------------------------------------ #

    def _get_ion_index(self, atoms_map: dict, target_symbol: str) -> str:
        """
        Given the symbol, returns the first matching ion index (str).
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
        Checks if elements have the same symbol, different indices,
        and have not been visited yet.
        """
        different_index = target_index != source_index
        same_symbol     = atoms_map.get(target_index) == symbol
        not_visited     = not visited_neighbors[target_index]
        return different_index and same_symbol and not_visited
