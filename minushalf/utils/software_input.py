"""
Utility to resolve the input files required for each supported software.
"""

from pathlib import Path


def get_input_filenames(
    software: str,
    input_file: str = None,
    indirect: bool = False,
) -> list:
    """
    Return the list of input files required for the specified software.

    Args:
        software: Name of the first-principles software.
        input_file: Needed input file for information retrieval.
        indirect: Whether the calculation requires the additional
            charge-density file.

    Returns:
        List of input filenames.
    """

    if software.upper() == "VASP":
        input_files = (
            ["INCAR", "POSCAR", "KPOINTS", "CHGCAR"]
            if indirect
            else ["INCAR", "POSCAR", "KPOINTS"]
        )

    elif software.upper() in ("QE", "QUANTUM_ESPRESSO"):

        pseudopotentials = _get_qe_pseudopotentials(input_file)

        input_files = [
            input_file,
            *pseudopotentials,
        ]

    else:
        raise ValueError(
            f"Unsupported software '{software}'. "
            "Supported software: VASP, QE."
        )

    return input_files


import re
from pathlib import Path

# Any of these starting a stripped line marks the end of ATOMIC_SPECIES.
_CARD_TERMINATOR_RE = re.compile(
    r"^(&\w+|ATOMIC_\w+|K_POINTS|CELL_PARAMETERS)\b", re.IGNORECASE
)


def _get_qe_pseudopotentials(input_file: str) -> list:
    """
    Read the pseudopotential filenames from the ATOMIC_SPECIES card
    of a Quantum ESPRESSO input file.

    Expected format:
        ATOMIC_SPECIES
        Si 28.086 Si.UPF
        C  12.011 C.UPF

    The pseudopotential filename is the third column.

    Args:
        input_file: Path to the Quantum ESPRESSO SCF input file.

    Returns:
        List of pseudopotential filenames.
    """
    scf_path = Path(input_file)
    if not scf_path.exists():
        raise FileNotFoundError(f"QE input file not found: {input_file}")

    lines = scf_path.read_text().splitlines()

    start = next(
        (i for i, line in enumerate(lines) if line.strip().upper() == "ATOMIC_SPECIES"),
        None,
    )
    if start is None:
        raise ValueError(f"Could not find ATOMIC_SPECIES card in {input_file}")
    
    pseudopotentials = []
    for line in lines[start + 1:]:
        stripped = line.strip()

        if not stripped:
            break  # blank line ends the card
        if stripped.startswith("!"):
            continue  # comment
        if _CARD_TERMINATOR_RE.match(stripped):
            break  # next card started

        columns = stripped.split()
        if len(columns) < 3:
            raise ValueError(
                f"Invalid ATOMIC_SPECIES entry in {input_file}: '{stripped}'"
            )
        pseudopotentials.append(columns[2])

    if not pseudopotentials:
        raise ValueError(
            f"No pseudopotentials found in ATOMIC_SPECIES card of {input_file}"
        )

    return pseudopotentials