"""
Utility to resolve output filenames for each software,
given the software name and optionally a software input file path.

For VASP: ignores input_name, returns hardcoded defaults.
For QE:   reads prefix and outdir from input_name, builds filenames.
"""
import re
import os
import glob


def get_output_filenames(software: str,
                         input_name: str = None,
                         base_path: str = None,
                         atom: str = None) -> dict:
    """
    Returns a dict of resolved output filenames for each factory method,
    given the software name and optionally a software input file path.

    Args:
        software       (str): software name, e.g. 'VASP' or 'QE'
        input_name (str): path to the software input file (required for QE)
        base_path      (str): base directory for relative paths
        atom           (str): chemical symbol e.g. 'Si'. If provided and
                              software is QE, resolves the UPF path for
                              that element from the ATOMIC_SPECIES card.

    Returns:
        filenames (dict): keys match factory method names, values are
                          resolved file paths.
            - "eigenvalues"
            - "fermi_energy"
            - "atoms_map"
            - "number_of_bands"
            - "number_of_kpoints"
            - "band_projection"
            - "nearest_neighbor"
            - "potential"
    """
    if software.upper() == "VASP":
        filenames = {
            "eigenvalues":       "EIGENVAL",
            "fermi_energy":      "vasprun.xml",
            "atoms_map":         "vasprun.xml",
            "number_of_bands":   "PROCAR",
            "number_of_kpoints": "PROCAR",
            "band_projection":   "PROCAR",
            "nearest_neighbor":  "OUTCAR",
            "potential":         "POTCAR",
        }
        if base_path:
            filenames = {
                k: os.path.join(base_path, v)
                for k, v in filenames.items()
            }

    elif software.upper() in ("QE", "QUANTUM_ESPRESSO"):
        if input_name is None:
            raise Exception(
                "QE requires --input-name to determine output file prefix."
            )
        prefix, outdir = _get_prefix_and_outdir(input_name)

        # Resolve outdir relative to the directory of the input file
        input_dir = os.path.dirname(os.path.abspath(input_name))
        if base_path:
            input_dir = base_path

        xml_file = os.path.join(outdir, f"{prefix}.xml")

        # Resolve UPF path for the requested atom, or None if not requested
        potential = _get_upf_for_atom(input_name, atom) if atom is not None else None

        filenames = {
            "eigenvalues":       xml_file,
            "fermi_energy":      xml_file,
            "atoms_map":         xml_file,
            "number_of_bands":   xml_file,
            "number_of_kpoints": xml_file,
            "band_projection":   _find_projwfc_up(input_dir, input_name),
            "nearest_neighbor":  xml_file,
            "potential":         potential,
        }
        if base_path:
            filenames = {
                k: os.path.join(base_path, v)
                for k, v in filenames.items()
            }

    else:
        raise Exception(
            f"input_name: unknown software '{software}'. "
            f"Supported: 'VASP', 'QE'."
        )

    return filenames


def _get_upf_for_atom(input_name_path: str, atom: str) -> str:
    """
    Reads a QE input file and returns the UPF filename for the
    requested chemical symbol, as declared in the ATOMIC_SPECIES card.

    Args:
        input_name_path (str): path to the QE input file (e.g. scf.in)
        atom            (str): chemical symbol to look up, e.g. 'Si'

    Returns:
        upf_path (str): path to the UPF file, resolved relative to the
                        directory of the input file.

    Raises:
        Exception: if the atom is not found in ATOMIC_SPECIES.
    """
    CARD_NAMES = {
        "ATOMIC_POSITIONS", "K_POINTS", "CELL_PARAMETERS",
        "CONSTRAINTS", "OCCUPATIONS", "ATOMIC_FORCES",
    }
    input_dir = os.path.dirname(os.path.abspath(input_name_path))

    with open(input_name_path) as fh:
        in_atomic_species = False
        for line in fh:
            stripped = line.strip()

            # Detect the start of the ATOMIC_SPECIES card
            if not in_atomic_species:
                if stripped.upper() == "ATOMIC_SPECIES":
                    in_atomic_species = True
                continue

            # Any card name or namelist start signals end of ATOMIC_SPECIES
            if stripped.startswith("&") or stripped.upper() in CARD_NAMES:
                break

            # Each line: Symbol  Mass  UPF_filename
            parts = stripped.split()
            if len(parts) >= 3 and parts[0] == atom:               
                return os.path.join(input_dir, parts[2])

    raise Exception(
        f"Atom '{atom}' not found in ATOMIC_SPECIES card of '{input_name_path}'."
    )


def _get_prefix_and_outdir(input_name_path: str) -> tuple:
    """
    Reads a QE input file and extracts 'prefix' and 'outdir'.

    Args:
        input_name_path (str): path to the QE input file (e.g. scf.in)

    Returns:
        (prefix, outdir) (tuple[str, str]):
            prefix — value of the prefix variable (e.g. 'AlN-wz')
            outdir — value of the outdir variable (e.g. './outdir/')
                     defaults to './' if not found
    """
    prefix_regex = re.compile(
        r"^\s*prefix\s*=\s*['\"]([^'\"]+)['\"]"
    )
    outdir_regex = re.compile(
        r"^\s*outdir\s*=\s*['\"]([^'\"]+)['\"]"
    )

    prefix = None
    outdir = "./"

    with open(input_name_path) as fh:
        for line in fh:
            if prefix is None:
                m = prefix_regex.match(line)
                if m:
                    prefix = m.group(1)

            outdir_match = outdir_regex.match(line)
            if outdir_match:
                outdir = outdir_match.group(1)

            if prefix is not None and outdir != "./":
                break

    if prefix is None:
        prefix = "pwscf"

    return prefix, outdir


def _find_projwfc_up(directory: str, input_name: str) -> Optional[str]:
    """
    Find the projwfc.x output file (*.projwfc_up).

    Lookup order:
        1. If the input file has a 'filproj = <prefix>' line, return
           '<directory>/<prefix>.projwfc_up'.
        2. Otherwise, search `directory` for a single *.projwfc_up file.
        3. If nothing is found, or several files match, return None.

    Args:
        directory (str): directory to search in
        input_name (str): name of the input file that may contain filproj

    Returns:
        path (str | None): path to the projwfc_up file, or None if it
        cannot be determined
    """

    _FILPROJ_RE = re.compile(r"\bfilproj\s*=\s*['\"]?([^'\",\s]+)", re.IGNORECASE)

    # 1. Look for filproj in the input file.
    input_path = os.path.join(directory, input_name)
    if os.path.isfile(input_path):
        with open(input_path, "r") as f:
            for line in f:
                line = line.split("!", 1)[0]  # drop comments
                match = _FILPROJ_RE.search(line)
                if match:
                    return os.path.join(directory, f"{match.group(1)}.projwfc_up")

    # 2. Fall back to searching the directory.
    matches = glob.glob(os.path.join(directory, "*.projwfc_up"))

    # 3. Exactly one match -> return it; zero or multiple -> None.
    if len(matches) == 1:
        return matches[0]
    return None