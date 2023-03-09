"""
Test constants
"""
import numpy as np
from minushalf.data.constants import Constants


def test_pi():
    """
    Test pi in constants module
    """
    constants = Constants()
    assert np.isclose(constants.pi_constant, 3.1415926548)


def test_trimming_exponent():
    """
    Test trimming exponent in constants module
    """
    constants = Constants()
    assert np.isclose(constants.trimming_exponent, 8)


def test_bohr_radius():
    """
    Test bohr radius in constants module
    """
    constants = Constants()
    assert np.isclose(constants.bohr_radius, 0.529177068)


def test_rydberg():
    """
    Test Rydberg constant in constants module
    """
    constants = Constants()
    assert np.isclose(constants.rydberg, 13.6058038)
