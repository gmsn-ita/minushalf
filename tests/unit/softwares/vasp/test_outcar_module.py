"""
Test Outcar module
The functions in this file obey the following rules:
test_(what the function is meant to do)_(chemical compound)
"""
import numpy as np
from minushalf.softwares.vasp import Outcar


def test_parse_distances_gan_3d(file_path):
    """
    Test if the distances are extracted correctly
    """
    filename = file_path("/gan-3d/OUTCAR")
    outcar = Outcar(filename)

    assert len(outcar.relative_distances['1']) == 4
    assert len(outcar.relative_distances['2']) == 4
    assert outcar.relative_distances['1'][0][0] == 2
    assert np.isclose(outcar.relative_distances['1'][0][1], 1.97)
    assert outcar.relative_distances['1'][-1][0] == 2
    assert np.isclose(outcar.relative_distances['1'][-1][1], 1.97)
    assert outcar.relative_distances['2'][0][0] == 1
    assert np.isclose(outcar.relative_distances['2'][0][1], 1.97)
    assert outcar.relative_distances['2'][-1][0] == 1
    assert np.isclose(outcar.relative_distances['2'][-1][1], 1.97)


def test_nearest_neighbor_gan_3d(file_path):
    """
    Test the nearest neighbor function
    """
    filename = file_path("/gan-3d/OUTCAR")
    outcar = Outcar(filename)

    assert np.isclose(outcar.nearest_neighbor_distance('1'), 1.97)
    assert np.isclose(outcar.nearest_neighbor_distance('2'), 1.97)


def test_parse_distances_sic_2d(file_path):
    """
    Test if the distances are extracted correctly
    """
    filename = file_path("/sic-2d/OUTCAR")
    outcar = Outcar(filename)

    assert len(outcar.relative_distances['1']) == 3
    assert len(outcar.relative_distances['2']) == 3
    assert outcar.relative_distances['1'][0][0] == 2
    assert np.isclose(outcar.relative_distances['1'][0][1], 1.79)
    assert outcar.relative_distances['1'][-1][0] == 2
    assert np.isclose(outcar.relative_distances['1'][-1][1], 1.78)
    assert outcar.relative_distances['2'][0][0] == 1
    assert np.isclose(outcar.relative_distances['2'][0][1], 1.79)
    assert outcar.relative_distances['2'][-1][0] == 1
    assert np.isclose(outcar.relative_distances['2'][-1][1], 1.79)


def test_nearest_neighbor_sic_2d(file_path):
    """
    Test the nearest neighbor function
    """
    filename = file_path("/sic-2d/OUTCAR")
    outcar = Outcar(filename)
    assert np.isclose(outcar.nearest_neighbor_distance('1'), 1.78)
    assert np.isclose(outcar.nearest_neighbor_distance('2'), 1.79)
