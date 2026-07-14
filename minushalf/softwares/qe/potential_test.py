"""
Test Potential class
"""
import numpy as np
from minushalf.softwares.qe.potential import Potential


def test_potential_al(file_path):
    """
    Test the potential parser for Aluminium
    """
    path = file_path('/Al/Al.upf')
    upf = Potential(path)

    assert upf.get_name() == "Al.upf" 
    assert upf.element == "Al"
    assert np.isclose(upf.z_valence, 3.0)
    assert np.isclose(upf.mesh_size, 1854)
    assert np.isclose(upf.r_grid[0], 0.0)
    assert np.isclose(upf.r_grid[-1], 18.53)
    assert np.isclose(len(upf.r_grid), 1854)
    assert np.isclose(upf.rab_grid[0], 0.01)
    assert np.isclose(upf.rab_grid[-1], 0.01)
    assert np.isclose(len(upf.rab_grid), 1854)
    assert np.isclose(upf.get_local_potential()[0], -4.5208790289)
    assert np.isclose(upf.get_local_potential()[-1], -4.8781370457)
    assert np.isclose(len(upf.get_local_potential()), 1854)

def test_potential_stringlist_al(file_path):
    """
    Test the potcar function
    to_stringlist for Aluminium
    """
    path = file_path('/Al/Al.upf')
    upf = Potential(path)

    with open(path, "r") as file:
        upf_generated_lines = upf.to_stringlist()
        for index, line in enumerate(file):
            assert line.strip() == upf_generated_lines[index].strip()



def test_potcar_c(file_path):
    """
    Test the potcar parser for Carbon
    """
    path = file_path('/C/C.upf')
    upf = Potential(path)

    assert upf.get_name() == "C.upf" 
    assert upf.element == "C"
    assert np.isclose(upf.z_valence, 4.0)
    assert np.isclose(upf.mesh_size, 1073)
    assert np.isclose(upf.r_grid[0], 0.0001519803275924194)
    assert np.isclose(upf.r_grid[-1], 100.3075063120137)
    assert np.isclose(len(upf.r_grid), 1073)
    assert np.isclose(upf.rab_grid[0], 1.899754094905242e-06)
    assert np.isclose(upf.rab_grid[-1], 1.253843828900171)
    assert np.isclose(len(upf.rab_grid), 1073)
    assert np.isclose(upf.get_local_potential()[0], -4.5208790289)
    assert np.isclose(upf.get_local_potential()[-1], -11.121142272)
    assert np.isclose(len(upf.get_local_potential()), 1073)


def test_potcar_stringlist_c(file_path):
    """
    Test the potcar function
    to_stringlist for Carbon
    """
    path = file_path('/C/C.upf')
    upf = Potential(path)

    with open(path, "r") as file:
        potential_generated_lines = upf.to_stringlist()
        for index, line in enumerate(file):
            assert line.strip() == potential_generated_lines[index].strip()

def test_potential_n(file_path):
    """
    Test the potential parser for Nitrogen
    """
    path = file_path('/N/N.upf')
    upf = Potential(path)

    assert upf.get_name() == "N.upf" 
    assert upf.element == "N"
    assert np.isclose(upf.z_valence, 5.0)
    assert np.isclose(upf.mesh_size, 1058)
    assert np.isclose(upf.r_grid[0], 0.0)
    assert np.isclose(upf.r_grid[-1], 10.57)
    assert np.isclose(len(upf.r_grid), 1058)
    assert np.isclose(upf.rab_grid[0], 0.01)
    assert np.isclose(upf.rab_grid[-1], 0.01)
    assert np.isclose(len(upf.rab_grid), 1058)
    assert np.isclose(upf.get_local_potential()[0], -17.397599866)
    assert np.isclose(upf.get_local_potential()[-1], -18.446066062)
    assert np.isclose(len(upf.get_local_potential()), 1058)

def test_potcar_stringlist_n(file_path):
    """
    Test the potcar function
    to_stringlist for Nitrogen
    """
    path = file_path('/N/N.upf')
    upf = Potential(path)

    with open(path, "r") as file:
        potential_generated_lines = upf.to_stringlist()
        for index, line in enumerate(file):
            assert line.strip() == potential_generated_lines[index].strip()

def test_potcar_lda_si(file_path):
    """
    Test the potcar using LDA for silicium
    """
    path = file_path('/Si/Si.upf')
    upf = Potential(path)

    assert upf.get_name() == "Si.upf" 
    assert upf.element == "Si"
    assert np.isclose(upf.z_valence, 4.0)
    assert np.isclose(upf.mesh_size, 431)
    assert np.isclose(upf.r_grid[0], 0.00130825992062)
    assert np.isclose(upf.r_grid[-1], 61.0041973233)
    assert np.isclose(len(upf.r_grid), 431)
    assert np.isclose(upf.rab_grid[0], 3.27064980156e-05)
    assert np.isclose(upf.rab_grid[-1], 1.52510493308)
    assert np.isclose(len(upf.rab_grid), 431)
    assert np.isclose(upf.get_local_potential()[0], -4.5208790289)
    assert np.isclose(upf.get_local_potential()[-1], -1.4440913966)
    assert np.isclose(len(upf.get_local_potential()), 431)


def test_potcar_stringlist_si_lda(file_path):
    """
    Test the potcar function
    to_stringlist for Silicium. The
    potcar uses LDA approximation.
    """
    path = file_path('/Si/Si.upf')
    upf = Potential(path)

    with open(path, "r") as file:
        potential_generated_lines = upf.to_stringlist()
        for index, line in enumerate(file):
            assert line.strip() == potential_generated_lines[index].strip()
