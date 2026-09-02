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
        if atom is not None:
            potential = _get_upf_for_atom(input_name, atom)
        else:
            potential = None

        filenames = {
            "eigenvalues":       xml_file,
            "fermi_energy":      xml_file,
            "atoms_map":         xml_file,
            "number_of_bands":   xml_file,
            "number_of_kpoints": xml_file,
            "band_projection":   _find_projwfc_up(input_dir),
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
    in_atomic_species = False

    with open(input_name_path) as fh:
        for line in fh:
            stripped = line.strip()

            # Detect the start of the ATOMIC_SPECIES card
            if stripped.upper() == "ATOMIC_SPECIES":
                in_atomic_species = True
                continue

            if not in_atomic_species:
                continue

            # Any card name or namelist start signals end of ATOMIC_SPECIES
            if stripped.startswith("&") or stripped.upper() in (
                "ATOMIC_POSITIONS", "K_POINTS",
                "CELL_PARAMETERS", "CONSTRAINTS",
                "OCCUPATIONS", "ATOMIC_FORCES"
            ):
                break

            # Each line: Symbol  Mass  UPF_filename
            parts = stripped.split()
            if len(parts) >= 3 and parts[0] == atom:
                upf_filename = parts[2]
                # Resolve relative to the input file's directory
                input_dir = os.path.dirname(os.path.abspath(input_name_path))
                return os.path.join(input_dir, upf_filename)

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


def _find_projwfc_up(directory: str) -> str:
    """
    Searches for a *.projwfc_up file in the given directory.

    This is a temporary solution until a dedicated projwfc.x input
    parser is implemented. It assumes there is exactly one projwfc_up
    file in the directory.

    Args:
        directory (str): directory to search in

    Returns:
        path (str): full path to the found .projwfc_up file

    Raises:
        Exception: if none or more than one file is found
    """
    matches = glob.glob(os.path.join(directory, "*.projwfc_up"))

    if not matches:
        raise Exception(
            f"input_name: could not find any *.projwfc_up file "
            f"in {directory}. Make sure projwfc.x has been run."
        )
    if len(matches) > 1:
        raise Exception(
            f"input_name: found multiple *.projwfc_up files in "
            f"{directory}: {matches}. Please remove the unused ones."
        )

    return matches[0]