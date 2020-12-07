"""
Parse valence orbital line in INP file.
"""


def parse_valence_orbitals(line: str) -> dict:
    """
    Parse valence orbital line in principal
    quantum number, angular momentum quantum number
    and eletronic occupation

        Args:
            line (str): line of imp file that represents
            a valence orbital
        Returns:
            A dictionary with fields n, l and eletronic occupation
        """
    parsed_line = line.split()
    orbital = {
        "n": int(parsed_line[0]),
        "l": int(parsed_line[1]),
        "occupation": [float(value) for value in parsed_line[2:]],
    }

    return orbital
