"""
Test ternary search
"""
import numpy as np
from minushalf.utils import ternary_search


def test_ternary_search_first():
    """
    Test ternary search with the cost function
    defined as follows:
        f(x) = -x**2 + 2x + 1
    """
    cost_function = lambda x: -x**2 + 2 * x + 1
    result = ternary_search(0, 10, cost_function)

    assert np.isclose(result, 1.0)


def test_ternary_search_second():
    """
    Test ternary search with the cost function
    defined as follows:
        f(x) = -3x**2 + 6x -2
    """
    cost_function = lambda x: -3 * x**2 + 6 * x - 2
    result = ternary_search(1.0, 10, cost_function)

    assert np.isclose(result, 1.0)


def test_ternary_search_third():
    """
    Test ternary search with the cost function
    defined as follows:
        f(x) = -2x**2 + 8x -10
    """
    cost_function = lambda x: -2 * x**2 + 8 * x - 10
    result = ternary_search(0.0, 2.0, cost_function)

    assert np.isclose(result, 2.0)


def test_ternary_search_fourth():
    """
    Test ternary search with the cost function
    defined as follows:
        f(x) = -4x**2 + 5x - 20
    """
    cost_function = lambda x: -4 * x**2 + 5 * x - 20
    result = ternary_search(0, 10, cost_function)

    assert np.isclose(result, 0.625)


def test_ternary_search_sixth():
    """
    Test ternary search with the cost function
    defined as follows:
        f(x) = -3x**2 + 4x -70
    """
    cost_function = lambda x: -3 * x**2 + 4 * x - 70
    result = ternary_search(0, 10, cost_function)

    assert np.isclose(result, 0.66666666666)
