"""
Reads procar file, an output of
VASP software
"""
import re
from itertools import islice
from loguru import logger
from minushalf.interfaces import BandProjectionFile


class Procar(BandProjectionFile):
    """
	Reads procar and store
	useful informations
	"""
    def __init__(self, filename: str):
        """
			Args:
				filename (str): name of the procar file in VASP
			Members:
				num_kpoints: Number of kpoints used in the simulation
				num_bands: Number of bands used in the simulation
				size_procar_header: Number of lines of the procar header in PROCAR file
		"""
        self.filename = filename
        self.num_kpoints, self.num_bands = self._get_num_kpoints_and_bands()
        self.size_procar_header = 3
        self.is_spin_polarized = self._is_spin_polarized()
        if self.is_spin_polarized:
            logger.error("The calculation is spin polarized")
            raise ValueError(
                "Minushalf CLI does not support calculations with spin polarized"
            )

    def _is_spin_polarized(self) -> bool:
        """
        Check if the calculation
        is spin polarized
        """
        header_regex = re.compile(
            r"^\s*#\s+of\s+k-points:\s+([0-9]+)\s+#\s+of\s+bands:\s+([0-9]+)")
        with open(self.filename, "r") as procar:

            for line in islice(procar, self.size_procar_header, None):
                header_match = header_regex.match(line)
                if header_match:
                    return True
        return False

    def _get_num_kpoints_and_bands(self) -> tuple:
        """
			Returns:
				(kpoints,band_number) (tuple): The number of kpoints
				used in the simulation and the number of bands used
				for each kpoint.
		"""
        regex = re.compile(
            r"^\s*#\s+of\s+k-points:\s+([0-9]+)\s+#\s+of\s+bands:\s+([0-9]+)")
        number_kpoints = None
        number_bands = None
        with open(self.filename) as procar:
            for line in procar:
                if regex.match(line):
                    number_kpoints = regex.match(line).group(1)
                    number_bands = regex.match(line).group(2)
                    return (int(number_kpoints), int(number_bands))

    def get_band_projection(self, kpoint: int, band_number: int):
        """
		Get the band projection for an specific kpoint and number of band
			Args:
				kpoint (int): Number of kpoints
				band_number (int): Number of the band
		"""

        projections = {}
        kpoint_regex = re.compile(r"^\s*k-point\s*([0-9]+)\s+")
        band_regex = re.compile(r"^\s*band\s*([0-9]+)")
        projections_regex = re.compile(
            r"^\s*([0-9]+)\s*([-+]?[0-9]*\.?[0-9]+)")

        with open(self.filename, "r") as procar:
            current_kpoint = None
            current_band = None

            for line in islice(procar, self.size_procar_header, None):

                kpoint_match = kpoint_regex.match(line)
                if kpoint_match:
                    current_kpoint = int(kpoint_match.group(1))

                band_match = band_regex.match(line)
                if band_match:
                    current_band = int(band_match.group(1))

                if current_kpoint == kpoint and current_band == band_number:
                    if projections_regex.match(line):
                        parsed_line = line.split()
                        projections["{}".format(parsed_line[0])] = [
                            float(elem) for elem in parsed_line[1:-1]
                        ]
                elif current_kpoint and current_kpoint > kpoint:
                    return projections

        return projections
