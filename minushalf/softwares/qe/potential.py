"""
Reads a UPF v2 pseudopotential file, an input file for
Quantum ESPRESSO software
"""
import os
import re
import numpy as np
from minushalf.softwares.potential_file import PotentialFile


class Potential(PotentialFile):
    """
    Parses a UPF v2 pseudopotential file and stores the local potential
    and radial mesh.
    """

    # Self-closing tag that holds all metadata as attributes:
    # <PP_HEADER ... mesh_size="1058" z_valence="5.00" element="N " .../>
    _HEADER_REGEX = re.compile(r"^\s*<PP_HEADER")

    # A data row: one or more scientific-notation or decimal numbers
    _DATA_ROW_REGEX = re.compile(
        r"^\s*[-+]?[0-9]*\.?[0-9]+(?:[EeDd][-+]?\d+)?"
    )

    # Attribute extraction helpers
    _ATTR_MESH_SIZE  = re.compile(r'mesh_size\s*=\s*"?\s*(\d+)"?')
    _ATTR_Z_VALENCE  = re.compile(r'z_valence\s*=\s*"?\s*([-+]?[0-9.E+\-]+)"?')
    _ATTR_ELEMENT    = re.compile(r'element\s*=\s*"([^"]+)"')

    def __init__(self, filename: str) -> None:
        """
        Args:
            filename (str): path to the UPF v2 file (e.g. 'N.upf')
        Members:
            filename    : stored path
            element     : chemical symbol parsed from PP_HEADER (str)
            z_valence   : valence charge from PP_HEADER (float)
            mesh_size   : number of radial grid points (int)
            r_grid      : radial grid r(i) in Bohr (np.ndarray, shape (mesh_size,))
            rab_grid    : integration weights dr(i) in Bohr (np.ndarray)
            potential   : local potential V_local(r) in Ry (np.ndarray)
        """
        self.filename = filename
        self.name = os.path.basename(self.filename)

        self.element, self.z_valence, self.mesh_size = (
            self._get_header_info()
        )
        self.r_grid  = self._get_block("PP_R")
        self.rab_grid = self._get_block("PP_RAB")
        self.potential = self._get_block("PP_LOCAL")

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
        """
        Reconstruct the UPF file as a list of strings, replacing only
        the PP_LOCAL block with the (possibly modified) self.potential.
        
        Returns:
            upf_lines (list[str]): full UPF file contents
        """
        local_open_re  = re.compile(r"^\s*<PP_LOCAL[^>]*>")
        local_close_re = re.compile(r"^\s*</PP_LOCAL\s*>")

        lines_out = []
        skip = False

        with open(self.filename, "r") as fh:
            for line in fh:
                if local_open_re.match(line):
                    # Write the original opening tag verbatim
                    lines_out.append(line)
                    # Inject the (modified) potential values
                    lines_out.extend(self._format_data_block(self.potential,
                                                             columns=4))
                    skip = True
                    continue
                if skip and local_close_re.match(line):
                    lines_out.append(line)
                    skip = False
                    continue
                if not skip:
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

    def _get_header_info(self) -> tuple[str, float, int]:
        """
        Parse the self-closing <PP_HEADER .../> tag to extract:
          - element symbol
          - z_valence
          - mesh_size

        The tag may span multiple lines, so lines are accumulated until
        the closing '/>' is found.

        Returns:
            (element, z_valence, mesh_size) (tuple[str, float, int])
        """
        header_lines = []
        inside = False

        with open(self.filename, "r") as fh:
            for line in fh:
                if not inside and self._HEADER_REGEX.match(line):
                    inside = True
                if inside:
                    header_lines.append(line)
                    if "/>" in line:
                        break

        if not header_lines:
            raise Exception(
                f"UPFFile parser could not find <PP_HEADER> in {self.filename}"
            )

        header_text = "".join(header_lines)

        element_match   = self._ATTR_ELEMENT.search(header_text)
        z_val_match     = self._ATTR_Z_VALENCE.search(header_text)
        mesh_size_match = self._ATTR_MESH_SIZE.search(header_text)

        if not (element_match and z_val_match and mesh_size_match):
            raise Exception(
                "UPFFile parser could not extract element, z_valence, or "
                f"mesh_size from PP_HEADER in {self.filename}"
            )

        return (
            element_match.group(1).strip(),
            float(z_val_match.group(1)),
            int(mesh_size_match.group(1)),
        )

    def _get_block(self, tag: str) -> np.ndarray:
        """
        Generic block reader: finds <tag ...> ... </tag> in the UPF file
        and returns all numeric values as a flat NumPy array.
        
        Args:
            tag (str): block name, e.g. 'PP_R', 'PP_RAB', 'PP_LOCAL'

        Returns:
            data (np.ndarray): flat float64 array of all values in the block
        """
        open_re  = re.compile(rf"^\s*<{re.escape(tag)}[^>]*>")
        close_re = re.compile(rf"^\s*</{re.escape(tag)}\s*>")

        values = []
        inside = False

        with open(self.filename, "r") as fh:
            for line in fh:
                if not inside and open_re.match(line):
                    inside = True
                elif inside and close_re.match(line):
                    break
                elif inside and self._DATA_ROW_REGEX.match(line):
                    values.extend(float(v) for v in line.split())

        if not values:
            raise Exception(
                f"UPFFile parser could not find block <{tag}> "
                f"in {self.filename}"
            )

        data = np.array(values, dtype=np.float64)

        if len(data) != self.mesh_size:
            raise Exception(
                f"UPFFile parser: block <{tag}> has {len(data)} values "
                f"but PP_HEADER declares mesh_size={self.mesh_size}"
            )

        return data

    @staticmethod
    def _format_data_block(array: np.ndarray, columns: int = 4) -> list[str]:
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

# ---------------------------------------------------------------------------
# Smoke-test
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    import sys

    filename = sys.argv[1] if len(sys.argv) > 1 else "/home/bruno-augusto/Desktop/QE_minushalf/QE/Al.upf"
    upf = Potential(filename)

    print(f"element    : {upf.element}")
    print(f"z_valence  : {upf.z_valence}")
    print(f"mesh_size  : {upf.mesh_size}")
    print(f"r_grid     : {upf.r_grid[:4]} ...")
    print(f"potential  : {upf.potential[:4]} ...")
    print(f"name       : {upf.name} ...")
    