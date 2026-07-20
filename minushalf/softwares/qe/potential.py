"""
Read a UPF v2 pseudopotential file
"""
import os
import numpy as np
import xml.etree.ElementTree as ET
from minushalf.softwares.potential_file import PotentialFile


class Potential(PotentialFile):
    """
    Parse a UPF v2 pseudopotential file and store the local potential
    and radial mesh.

    UPF v2 file structure (relevant blocks):
    -----------------------------------------
    <PP_HEADER ... mesh_size="1058" z_valence="5.00" element="N" .../>

    <PP_MESH>
      <PP_R type="real" size="1058" columns="8">
        0.0000    0.0100    0.0200 ...      ← radial grid r(i) in Bohr
      </PP_R>
      <PP_RAB type="real" size="1058" columns="8">
        ...                                 ← integration weights dr(i)
      </PP_RAB>
    </PP_MESH>

    <PP_LOCAL type="real" size="1058" columns="4">
      -1.7028214881E+01  -1.7027517148E+01  ...   ← V_local(r) in Ry
    </PP_LOCAL>

    Values in PP_LOCAL are V_local(r) in Ry (QE convention), sampled on
    the radial grid. The radial grid PP_R and weights
    PP_RAB are stored alongside so that callers can perform a Fourier
    transform if needed.
    """

    def __init__(self, filename: str) -> None:
        """
        Args:
            filename (str): path to the UPF v2 file (e.g. 'N.upf')
        Members:
            filename    : stored path
            name        : name of the upf file
            element     : chemical symbol parsed from PP_HEADER (str)
            z_valence   : valence charge from PP_HEADER (float)
            mesh_size   : number of radial grid points (int)
            r_grid      : radial grid r(i) in Bohr (np.ndarray, shape (mesh_size,))
            rab_grid    : integration weights dr(i) in Bohr (np.ndarray)
            potential   : local potential V_local(r) in Ry (np.ndarray)
        """
        self.filename = filename
        self._root = ET.parse(filename).getroot()
        self.name = os.path.basename(self.filename)

        self.element, self.z_valence, self.mesh_size = (
            self._get_header_info()
        )
        self.r_grid   = self._get_block("PP_R")
        self.rab_grid = self._get_block("PP_RAB")
        self.potential = self._get_block("PP_LOCAL")

    # ------------------------------------------------------------------ #
    #  Public interface — mirrors Potcar                                   #
    # ------------------------------------------------------------------ #

    def get_local_potential(self) -> np.ndarray:
        """
        Returns the local potential array V_local(r) sampled on the
        radial grid.
 
        Returns:
            potential (np.ndarray): V_local(r) values in Ry,
                                    shape (mesh_size,)
        """
        return self.potential

    def get_name(self) -> str:
        """
        Returns potential file name tag.
        """
        return self.name

    def get_maximum_module_wave_vector(self) -> None:
        """
        Not applicable for QE: k_max is not used by the QE correction
        workflow and is not provided in the UPF file.

        If an estimate is ever needed, the Nyquist limit from the radial
        grid spacing can be used:
            dr_min = min(diff(r_grid[r_grid > 0]))
            k_max  ≈ π / dr_min   [Bohr⁻¹]
        """
        pass

    def to_stringlist(self) -> list:
        lines_out = []
        skip = False

        with open(self.filename, "r") as fh:
            for line in fh:
                stripped = line.strip()
                if stripped.startswith("<PP_LOCAL") and stripped.endswith(">"):
                    lines_out.append(line)
                    lines_out.extend(self._format_data_block(self.potential,
                                                            columns=4))
                    skip = True
                elif skip and stripped == "</PP_LOCAL>":
                    lines_out.append(line)
                    skip = False
                elif not skip:
                    lines_out.append(line)

        return lines_out

    def to_file(self, filename: str) -> None:
        """
        Write the (possibly modified) UPF file to disk.

        Args:
            filename (str): output file path
        """
        lines = self.to_stringlist()
        with open(filename, "w") as fh:
            fh.writelines(lines)

    # ------------------------------------------------------------------ #
    #  Private helpers                                                     #
    # ------------------------------------------------------------------ #

    def _get_header_info(self) -> tuple:
        """
        Extract element, z_valence and mesh_size from the
        <PP_HEADER .../> attributes using ElementTree.

        Returns:
            (element, z_valence, mesh_size) (tuple[str, float, int])
        """
        header = self._root.find(".//PP_HEADER")
        if header is None:
            raise Exception(
                f"UPFFile parser could not find <PP_HEADER> in {self.filename}"
            )

        element   = header.get("element")
        z_valence = header.get("z_valence")
        mesh_size = header.get("mesh_size")

        if not (element and z_valence and mesh_size):
            raise Exception(
                "UPFFile parser could not extract element, z_valence, or "
                f"mesh_size from PP_HEADER in {self.filename}"
            )

        return (
            element.strip(),
            float(z_valence),
            int(mesh_size),
        )

    def _get_block(self, tag: str) -> np.ndarray:
        """
        Generic block reader: finds the element with the given tag anywhere
        in the tree and returns all its numeric text content as a flat
        NumPy array.

        Args:
            tag (str): block name, e.g. 'PP_R', 'PP_RAB', 'PP_LOCAL'

        Returns:
            data (np.ndarray): flat float64 array of all values in the block
        """
        element = self._root.find(f".//{tag}")
        if element is None or not element.text:
            raise Exception(
                f"UPFFile parser could not find block <{tag}> "
                f"in {self.filename}"
            )

        data = np.array(element.text.split(), dtype=np.float64)

        if len(data) != self.mesh_size:
            raise Exception(
                f"UPFFile parser: block <{tag}> has {len(data)} values "
                f"but PP_HEADER declares mesh_size={self.mesh_size}"
            )

        return data

    @staticmethod
    def _format_data_block(array: np.ndarray, columns: int = 4) -> list:
        """
        Format a NumPy array back into UPF-style lines of `columns`
        values each, using scientific notation matching QE's output.

        Args:
            array   (np.ndarray): values to format
            columns (int):        values per line (default 4, as in PP_LOCAL)

        Returns:
            lines (list[str]): formatted lines including newlines
        """
        lines = []
        for i in range(0, len(array), columns):
            chunk = array[i:i + columns]
            line = "   ".join(f"{v:18.10E}" for v in chunk)
            lines.append(f"  {line}\n")
        return lines
    