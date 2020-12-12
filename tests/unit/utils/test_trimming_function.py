"""
Test trimming function
"""
import numpy as np
from minushalf.utils import trimming_function


def test_trimming_first():
    """
    Test the trimming function with the following parameters:
        cut: 3.0
        amplitude:1.0
        radius: 3.0
        ion_potential: 0.77e-5
        atom_potential: 0.64e-4
    """
    res = trimming_function(3.0, 0.77e-5, 0.64e-4, 3.0, 1.0)
    assert np.isclose(res, 0)


def test_trimming_second():
    """
    Test the trimming function with the following parameters:
        cut: 3.0
        amplitude:1.0
        radius: 2.9
        ion_potential: 0.77e-5
        atom_potential: 0.64e-4
    """
    res = trimming_function(2.9, 0.77e-5, 0.64e-4, 3.0, 1.0)
    assert np.isclose(res, -1.911990822901309e-05)


def test_trimming_third():
    """
    Test the trimming function with the following parameters:
        cut: 3.0
        amplitude:2.0
        radius: 2.9
        ion_potential: 0.77e-5
        atom_potential: 0.64e-4
    """
    res = trimming_function(2.9, 0.77e-5, 0.64e-4, 3.0, 2.0)
    assert np.isclose(res, -3.823981645802618e-05)


def test_trimming_fourth():
    """
    Test the trimming function with the following parameters:
        cut: 1.55
        amplitude:1.54
        radius: 0.23e-2
        ion_potential: 1.2
        atom_potential: 0.34
    """
    res = trimming_function(0.23e-2, 1.2, 0.34, 1.55, 1.54)
    assert np.isclose(res, 33.55492614903707)


def test_trimming_fifth():
    """
    Test the trimming function with the following parameters:
        cut: 1.55
        amplitude:1.54
        radius: 100
        ion_potential: 1.2
        atom_potential: 0.34
    """
    res = trimming_function(100, 1.2, 0.34, 1.55, 1.54)
    assert np.isclose(res, 0)


def test_trimming_sixth():
    """
    Test the trimming function with the following parameters:
        cut: 1.55
        amplitude:1.54
        radius: 1
        ion_potential: 1.2e-23
        atom_potential: 0.34e-34
    """
    res = trimming_function(1, 1.2e-23, 0.34e-34, 1.55, 1.54)
    assert np.isclose(res, 4.273004874113689e-22)


def test_trimming_seventh():
    """
    Test the trimming function with the following parameters:
        cut: 3.76
        amplitude:0.34
        radius: 2.95
        ion_potential: 1.643
        atom_potential: 1.643
    """
    res = trimming_function(2.95, 1.643, 1.643, 3.76, 0.34)
    assert np.isclose(res, 0)


def test_trimming_eighth():
    """
    Test the trimming function with the following parameters:
        cut: 3.76
        amplitude:1.0
        radius: 2.95
        ion_potential: 1.643
        atom_potential: 1.835
    """
    res = trimming_function(2.95, 1.643, 1.835, 3.76, 1.0)
    assert np.isclose(res, -3.0556886231671694)


def test_trimming_nineth():
    """
    Test the trimming function with the following parameters:
        cut: 3.76
        amplitude:1.0
        radius: 0
        ion_potential: 1.643
        atom_potential: 1.835
    """
    res = trimming_function(0, 1.643, 1.835, 3.76, 1.0)
    assert np.isclose(res, -4.8645015256834165)


def test_trimming_tenth():
    """
    Test the trimming function with the following parameters:
        cut: 4.55
        amplitude:1.0
        radius: 4.32
        ion_potential: 12e-23
        atom_potential: 15e-34
    """
    res = trimming_function(4.32, 12e-23, 15e-34, 4.55, 1.0)
    assert np.isclose(res, 1.1912043652717141e-22)
