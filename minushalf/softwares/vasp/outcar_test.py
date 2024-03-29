"""
Test Outcar module
The functions in this file obey the following rules:
test_(what the function is meant to do)_(chemical compound)
"""
import pytest   
import numpy as np
from minushalf.softwares.vasp.outcar import Outcar


def test_parse_distances_gan_3d(file_path):
    """
    Test if the distances are extracted correctly
    """
    filename = file_path("/gan-3d/OUTCAR")
    outcar = Outcar(filename)

    assert len(outcar.relative_distances["1"]) == 4
    assert len(outcar.relative_distances["2"]) == 4
    assert outcar.relative_distances["1"][0][0] == 2
    assert np.isclose(outcar.relative_distances["1"][0][1], 1.97)
    assert outcar.relative_distances["1"][-1][0] == 2
    assert np.isclose(outcar.relative_distances["1"][-1][1], 1.97)
    assert outcar.relative_distances["2"][0][0] == 1
    assert np.isclose(outcar.relative_distances["2"][0][1], 1.97)
    assert outcar.relative_distances["2"][-1][0] == 1
    assert np.isclose(outcar.relative_distances["2"][-1][1], 1.97)


def test_nearest_neighbor_gan_3d(file_path):
    """
    Test the nearest neighbor function
    """
    filename = file_path("/gan-3d/OUTCAR")
    outcar = Outcar(filename)

    assert np.isclose(outcar.nearest_neighbor_distance("1"), 1.97)
    assert np.isclose(outcar.nearest_neighbor_distance("2"), 1.97)


def test_parse_distances_sic_2d(file_path):
    """
    Test if the distances are extracted correctly
    """
    filename = file_path("/sic-2d/OUTCAR")
    outcar = Outcar(filename)

    assert len(outcar.relative_distances["1"]) == 3
    assert len(outcar.relative_distances["2"]) == 3
    assert outcar.relative_distances["1"][0][0] == 2
    assert np.isclose(outcar.relative_distances["1"][0][1], 1.79)
    assert outcar.relative_distances["1"][-1][0] == 2
    assert np.isclose(outcar.relative_distances["1"][-1][1], 1.78)
    assert outcar.relative_distances["2"][0][0] == 1
    assert np.isclose(outcar.relative_distances["2"][0][1], 1.79)
    assert outcar.relative_distances["2"][-1][0] == 1
    assert np.isclose(outcar.relative_distances["2"][-1][1], 1.79)


def test_nearest_neighbor_sic_2d(file_path):
    """
    Test the nearest neighbor function
    """
    filename = file_path("/sic-2d/OUTCAR")
    outcar = Outcar(filename)
    assert np.isclose(outcar.nearest_neighbor_distance("1"), 1.78)
    assert np.isclose(outcar.nearest_neighbor_distance("2"), 1.79)


def test_number_equal_neighbors_sic_2d(file_path):
    """
    Test the nearest equal neighbor function
    """
    filename = file_path("/sic-2d/OUTCAR")
    fake_atoms_map = {"1": "Si", "2": "Si"}
    real_atoms_map = {"1": "Si", "2": "C"}
    outcar = Outcar(filename)
    assert outcar.number_of_equal_neighbors(fake_atoms_map, "Si") == 1
    assert outcar.number_of_equal_neighbors(real_atoms_map, "Si") == 0

def test_nearest_neighbor_aln_3d(file_path):
    """
    Test the nearest neighbor funcion
    """
    filename = file_path("/aln-3d/OUTCAR")
    outcar = Outcar(filename)
    assert np.isclose(outcar.nearest_neighbor_distance("4"), 1.90)

def test_nearest_neighbor_rus2(file_path):
    """
    Test the nearest neighbor funcion
    """
    filename = file_path("/rus2/OUTCAR")
    outcar = Outcar(filename)
    assert np.isclose(outcar.nearest_neighbor_distance("10"), 2.21)

@pytest.mark.xfail
def test_nearest_neighbor_rus2_fail(file_path):
    """
    Test the nearest neighbor funcion
    """
    filename = file_path("/rus2/OUTCAR")
    outcar = Outcar(filename)
    assert np.isclose(outcar.nearest_neighbor_distance("0"), 2.21)

def test_nearest_neighbor_pbte(file_path):
    """
    Test the nearest neighbor funcion
    """
    filename = file_path("/pbte/OUTCAR")
    outcar = Outcar(filename)
    assert np.isclose(outcar.nearest_neighbor_distance("2"), 3.28)

