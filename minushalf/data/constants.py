"""
Physical constants
"""


class Constants:
    """
    Class for physical constants used in the program. Contains:

        pi: About 3.1415

        trimming exponent: exponent used in the trimming function

        bohr_radius: The Bohr radius is a physical constant, equal to
                     the most probable distance between the nucleus and
                     the electron in a hydrogen atom in its ground state.

        rydberg: In spectroscopy, the Rydberg constant, symbol for heavy
                 atoms or for hydrogen, named after the Swedish physicist
                 Johannes Rydberg, is a physical constant relating to the
                 electromagnetic spectra of an atom
    """
    def __init__(self) -> None:
        self.pi_constant = 3.1415926548
        self.trimming_exponent = 8
        self.bohr_radius = 0.529177068
        self.rydberg = 13.6058038
