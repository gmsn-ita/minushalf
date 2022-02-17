"""
Reads Outcar file, an output of
VASP software
"""
import re
from collections import defaultdict


class Outcar():
    """
    Reads output informations stored
    in OUTCAR file
    """
    def __init__(self, filename: str):
        """
            Args:
                filename (str): name of the OUTCAR file in VASP

            Members:
                relative_distances (defaultdict(list)): An dictionary where
                the keys are the ion index given by VASP and the values are lists of tuples
                containing the index other ion and the relative distance
                respectively.
        """
        self.filename = filename
        self.relative_distances = self._get_distances()

    def nearest_neighbor_distance(self, ion_index: str) -> float:
        """
        Given the ion index, it returns the distance of the nearest neighbor
        to this ion

            Args:
                ion_index (str): The index of the ion given by VASP.

            Returns:
                nearest_neighbor_distance (float): The distance of the nearest neighbor.
        """
        distances = [dist for _, dist in self.relative_distances[ion_index]]
        return min(distances)

    def _get_ion_index(self, atoms_map: dict, target_symbol: str) -> int:
        """
        Given the symbol, returns the ion index
        """
        for index, symbol in atoms_map.items():
            if target_symbol == symbol:
                return index

    def _is_neighbor_equal(
        self,
        target_index: int,
        source_index: int,
        visited_neighbors: dict,
        atoms_map: dict,
        symbol: str,
    ) -> bool:
        """
        checks if elements have the same symbol and different indexes
        """
        different_index = target_index != source_index
        same_symbol = atoms_map[target_index] == symbol
        not_visited = not visited_neighbors[target_index]

        return different_index and same_symbol and not_visited

    def number_of_equal_neighbors(self, atoms_map: dict, symbol: str) -> int:
        """
        Given an map that links atoms symbols with it's index
        this function returns the number of neighbors of the atom with
        equal symbol but different indexes.

            Args:
                atoms_map (dict): Map the atoms index to their symbol.
                symbom (str): The symbol of the target atom.

            Returns:
                number_equal_neighbors (int): Returns the number of neighbors with
                                        same symbol but different indexes.
        """
        ion_index = self._get_ion_index(atoms_map, symbol)
        number_equal_neighbors = 0
        visited_neighbors = defaultdict(bool)  ## Already listed neighbors

        for index, _ in self.relative_distances[ion_index]:
            index = str(index)
            if self._is_neighbor_equal(index, ion_index,
                                       visited_neighbors, atoms_map, symbol):
                visited_neighbors[index] = True
                number_equal_neighbors += 1

        return number_equal_neighbors

    def _get_distances(self) -> defaultdict(list):
        """
            Returns:
                relative_distances (defaultdict(list)): An dictionary where
                the keys are the ion index given by VASP and the values are lists of tuples
                containing the index other ion and the relative distance
                respectively.
        """
        relative_distances = defaultdict(list)
        start_capture_regex = re.compile(
            r"\s*ion\s+position\s+nearest\s+neighbor\s+table")
        distances_line_regex = re.compile(
            r"\s+([0-9]).*(?=-)-\s+([0-9]\s+[0-9]*\.[0-9]+\s*)+")
        detect_break_line = re.compile(r"\s*[0-9]\s+[0-9]*\.[0-9]+\s*")

        with open(self.filename) as outcar:

            ## Skip initial lines
            for line in outcar:
                if start_capture_regex.match(line):
                    break
            ## Read informations
            ion_index = ""
            for line in outcar:

                if distances_line_regex.match(line) or detect_break_line.match(line):
                    if distances_line_regex.match(line):
                        ion_index = distances_line_regex.match(line).group(1)
                     
                    distance_line = ""
                    if '-' in line:
                        distance_line = line.split("-")[1]
                    else:
                        distance_line = line

                    ion_relative_distances_and_index = re.findall(
                            r"\s*[0-9]\s+[0-9]*\.[0-9]+\s*",
                        distance_line)

                    for element in ion_relative_distances_and_index:
                        element = element.rstrip()
                        index, distance = int(element.split()[0]), float(
                            element.split()[1])
                        relative_distances[ion_index].append((index, distance))
                else:
                    break

        return relative_distances
