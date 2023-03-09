"""
Test parse cut function
"""
import numpy as np
from minushalf.utils.parse_cut import parse_cut


def test_cut_range():
    """
    Test parse_cut with range
    """
    cut_range = "3.51:0.01:3.80"
    cuts = parse_cut(cut_range)
    expected_cuts = np.arange(3.51, 3.80, 0.01)
    for index, cut in enumerate(expected_cuts):
        assert np.isclose(cuts[index], cut)


def test_cut_range():
    """
    Test parse_cut with range
    """
    cut_str = "3.51"
    cut = parse_cut(cut_str)
    assert np.isclose(cut, 3.51)
