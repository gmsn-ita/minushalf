"""
Test initial cut guess module
"""
import numpy as np
from minushalf.utils import CutInitialGuess


def test_3d_method():
    """
    Test the method of guess for 3d compunds
    """
    cut_guess = CutInitialGuess()
    guess_1 = cut_guess.guess(0, "3d")
    assert np.isclose(guess_1, 0.15)
    guess_1 = cut_guess.guess(1, "3d")
    assert np.isclose(guess_1, 0.84 * 1.88973 + 0.15)
    guess_1 = cut_guess.guess(2, "3d")
    assert np.isclose(guess_1, 0.84 * 1.88973 * 2 + 0.15)
