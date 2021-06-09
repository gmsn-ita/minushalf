"""
Test correct potential fourier transform
"""
import pytest
import numpy as np
from minushalf.utils import correct_potential_fourier_transform


def test_correct_potential_fourier_transform_first():
    """
    Test correct_potential_fourier_transform with the
    following parameters:

    coefficient: 0
    k: 0
    rays: np.arange(0,5,0.1)
    occupation_potential: np.arange(20,25,0.1)
    cut: 3.0

    """
    res = correct_potential_fourier_transform(np.array([0]),
                                              np.array([0], dtype=float),
                                              np.arange(0, 5, 0.1),
                                              np.arange(20, 25, 0.1), 3.0)
    assert np.isclose(res[0], 92.23450000000014)


def test_correct_potential_fourier_transform_second():
    """
    Test correct_potential_fourier_transform with the
    following parameters:

    coefficient: 1
    k: 0
    rays: np.arange(0,5,0.1)
    occupation_potential: np.arange(20,25,0.1)
    cut: 3.0

    """
    res = correct_potential_fourier_transform(np.array([1]),
                                              np.array([0], dtype=float),
                                              np.arange(0, 5, 0.1),
                                              np.arange(20, 25, 0.1), 3.0)
    assert np.isclose(res[0], 93.23450000000014)


def test_correct_potential_fourier_transform_third():
    """
    Test correct_potential_fourier_transform with the
    following parameters:

    coefficient: 1
    k: 2e-1
    rays: np.arange(0,5,0.1)
    occupation_potential: np.arange(20,25,0.1)
    cut: 3.0

    """
    res = correct_potential_fourier_transform(np.array([1]), np.array([2e-1]),
                                              np.arange(0, 5, 0.1),
                                              np.arange(20, 25, 0.1), 3.0)
    assert np.isclose(res[0], 92.49911868134099)


def test_correct_potential_fourier_transform_fourth():
    """
    Test correct_potential_fourier_transform with the
    following parameters:

    coefficient: 4
    k: 3
    rays: np.arange(12,19,0.01)
    occupation_potential: np.arange(70,79,0.01)
    cut: 13.6

    """
    res = correct_potential_fourier_transform(np.array([4]), np.array([3]),
                                              np.arange(12, 19, 0.01),
                                              np.arange(70, 79, 0.01), 13.6)
    assert np.isclose(res[0], 110.19498579254733)


def test_correct_potential_fourier_transform_fifth():
    """
    Test correct_potential_fourier_transform with the
    following parameters:

    coefficient: 0
    k: 3
    rays: np.arange(12,19,0.01)
    occupation_potential: np.arange(70,79,0.01)
    cut: 13.6

    """
    res = correct_potential_fourier_transform(np.array([0]), np.array([3]),
                                              np.arange(12, 19, 0.01),
                                              np.arange(70, 79, 0.01), 13.6)
    assert np.isclose(res[0], 106.19498579254733)


@pytest.mark.xfail
def test_correct_potential_fourier_transform_sixth():
    """
    Test correct_potential_fourier_transform with the
    following parameters:

    coefficient: 0
    k: 3
    rays: np.arange(12,19,0.01)
    occupation_potential: np.arange(70,79,0.01)
    cut: 11.6

    """
    correct_potential_fourier_transform(np.array([0]), np.array([3]),
                                        np.arange(12, 19, 0.01),
                                        np.arange(70, 79, 0.01), 11.6)


def test_correct_potential_fourier_transform_seventh():
    """
    Test correct_potential_fourier_transform with the
    following parameters:

    coefficient: 0
    k: 67
    rays: np.arange(0,12,0.01)
    occupation_potential: np.arange(70,79,0.01)
    cut: 1.23

    """
    res = correct_potential_fourier_transform(np.array([0]), np.array([67]),
                                              np.arange(0, 12, 0.01),
                                              np.arange(70, 79, 0.01), 1.23)
    assert np.isclose(res[0], 0.013212154569150683)


def test_correct_potential_fourier_transform_eighth():
    """
    Test correct_potential_fourier_transform with the
    following parameters:

    coefficient: 0
    k: 67
    rays: np.arange(0,12,0.1)
    occupation_potential: np.arange(70,79,0.1)
    cut: 1.23

    """
    res = correct_potential_fourier_transform(np.array([0]), np.array([67]),
                                              np.arange(0, 12, 0.1),
                                              np.arange(70, 79, 0.1), 1.23)
    assert np.isclose(res[0], -0.017530596213630182)


def test_correct_potential_fourier_transform_nineth():
    """
    Test correct_potential_fourier_transform with the
    following parameters:

    coefficient: 0
    k: 67
    rays: np.arange(0,12,0.2)
    occupation_potential: np.arange(70,79,0.2)
    cut: 1.23

    """
    res = correct_potential_fourier_transform(np.array([0]), np.array([67]),
                                              np.arange(0, 12, 0.2),
                                              np.arange(70, 79, 0.2), 1.23)
    assert np.isclose(res[0], 0.3972120670340923)


def test_correct_potential_fourier_transform_tenth():
    """
    Test correct_potential_fourier_transform with the
    following parameters:

    coefficient: 0
    k: 67
    rays: np.arange(0,12,0.4)
    occupation_potential: np.arange(70,79,0.4)
    cut: 1.23

    """
    res = correct_potential_fourier_transform(np.array([0]), np.array([67]),
                                              np.arange(0, 12, 0.4),
                                              np.arange(70, 79, 0.4), 1.23)
    assert np.isclose(res[0], 0.32399512659710605)
