"""
Reads Outcar file, an output of
VASP software
"""
import re
import math
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
        distances = self.relative_distances[ion_index]
        nearest_distance = math.inf
        for distance in distances:
            nearest_distance = min(distance[1], nearest_distance)
        return nearest_distance

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
        ion_index = None
        for key, value in atoms_map.items():
            if value == symbol:
                ion_index = key
                break

        number_equal_neighbors = 0
        ## Already listed neighbors
        visited_neighbors = {key: False for key in atoms_map.keys()}
        for element in self.relative_distances[ion_index]:
            if element[0] != ion_index and atoms_map[str(
                    element[0])] == symbol and not visited_neighbors[str(
                        element[0])]:
                visited_neighbors[str(element[0])] = True
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
        with open(self.filename) as outcar:
            start_capture = False
            for line in outcar:

                if start_capture and distances_line_regex.match(line):
                    ion_index = distances_line_regex.match(line).group(1)

                    ion_relative_distances_line = line.split("-")[1]
                    ion_relative_distances_and_index = re.findall(
                        r"[0-9]\s+[0-9]*\.[0-9]+\s*",
                        ion_relative_distances_line)

                    for element in ion_relative_distances_and_index:
                        element_cleaned = element.rstrip()
                        index = int(element_cleaned.split()[0])
                        distance = float(element_cleaned.split()[1])
                        relative_distances[ion_index].append((index, distance))
                elif start_capture:
                    break

                if start_capture_regex.match(line):
                    start_capture = True

        return relative_distances
