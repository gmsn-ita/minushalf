"""
List atomic orbitals and
their respective groups
"""
from enum import Enum, unique


@unique
class OrbitalType(Enum):
    """
    Enum type for orbital type. Indices are basically the azimuthal quantum
    number, l.
    """

    s = 0
    p = 1
    d = 2
    f = 3

    def __str__(self):
        return str(self.name)


@unique
class Orbital(Enum):
    """
    Enum type for specific orbitals. The indices are basically the order in
    which the orbitals are reported in VASP and has no special meaning.
    """

    s = 0
    py = 1
    pz = 2
    px = 3
    dxy = 4
    dyz = 5
    dz2 = 6
    dxz = 7
    dx2 = 8
    f_3 = 9
    f_2 = 10
    f_1 = 11
    f0 = 12
    f1 = 13
    f2 = 14
    f3 = 15

    def __str__(self):
        return str(self.name)
