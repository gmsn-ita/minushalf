"""
Test orbital script
in utils module
"""
from minushalf.utils import ORBITAL_GROUP, ORBITALS


def test_orbitals():
    """
    Verify if the list returned
    contains all the orbitals
    """
    expect_orbitals = ["s", "px", "py", "pz",
                       "dxy", "dyz", "dz2", "dxz", "dx2y2"]

    assert expect_orbitals.__eq__(ORBITALS)


def test_orbital_group():
    """
    Verify if the list returned
    contains all the groups for each orbital
    """
    expect_orbital_group = set(["s", "p", "d"])

    assert expect_orbital_group.__eq__(ORBITAL_GROUP)
