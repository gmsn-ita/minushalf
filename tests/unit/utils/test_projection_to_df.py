"""
Test projection to df method
"""
import pytest
import numpy as np
from minushalf.utils import projection_to_df


def test_projection_to_df_first():
    """
    Test the method with the following input
    {"In": [1, 2, 3, 4, 5], "Sb": [1, 2, 3, 4, 5, 6]}
    """
    projections = {"In": [1, 2, 3, 4, 5], "Sb": [1, 2, 3, 4, 5, 6]}

    data_frame = projection_to_df(projections)
    assert np.isclose(data_frame["d"]["In"], 14.0)
    assert np.isclose(data_frame["d"]["Sb"], 31.0)
    assert np.isclose(data_frame["p"]["In"], 25.0)
    assert np.isclose(data_frame["p"]["Sb"], 25.0)
    assert np.isclose(data_frame["s"]["In"], 3.0)
    assert np.isclose(data_frame["s"]["Sb"], 3.0)


def test_projection_to_df_second():
    """
    Test the method with the following input
    {
        "In": [1, 2, 3, 4, 5],
        "Sb": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    }
    """
    projections = {
        "In": [1, 2, 3, 4, 5],
        "Sb": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    }

    data_frame = projection_to_df(projections)
    assert np.isclose(data_frame["d"]["In"], 7.0)
    assert np.isclose(data_frame["d"]["Sb"], 50.0)
    assert np.isclose(data_frame["p"]["In"], 13.0)
    assert np.isclose(data_frame["p"]["Sb"], 13.0)
    assert np.isclose(data_frame["s"]["In"], 1.0)
    assert np.isclose(data_frame["s"]["Sb"], 1.0)
    assert np.isclose(data_frame["f"]["In"], 0.0)
    assert np.isclose(data_frame["f"]["Sb"], 14.0)


@pytest.mark.xfail
def test_projection_to_df_third():
    """
    Test the method with the following input
    {
        "In": [],
        "Sb": []
    }
    """
    projections = {"In": [], "Sb": []}

    projection_to_df(projections)


@pytest.mark.xfail
def test_projection_to_df_fourth():
    """
    Test the method with the following input
    {
        "In": [1, 2, 3, 4, 5],
        "Sb": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10,11,12,13,14,15,16,17,18]
    }
    """
    projections = {
        "In": [1, 2, 3, 4, 5],
        "Sb": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18]
    }

    projection_to_df(projections)
