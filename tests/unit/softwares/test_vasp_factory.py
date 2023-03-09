"""
Test vasp factory module
"""
from yaml import compose_all
from minushalf import commands
import numpy as np
from minushalf.softwares.vasp.procar import Procar
from minushalf.softwares.vasp.potcar import Potcar
from minushalf.softwares.vasp.runner import VaspRunner
from minushalf.softwares.vasp_factory import Vasp


def test_get_atoms_map(file_path):
    """
    Test get atoms map function
    """
    base_path = file_path("/gec-2d/")
    factory = Vasp()
    atoms_map = factory.get_atoms_map(base_path=base_path)
    assert atoms_map["1"] == "Ge"
    assert atoms_map["2"] == "C"


def test_get_band_projection_class(file_path):
    """
    Test get band projection class function
    """
    base_path = file_path("/gec-2d/")
    factory = Vasp()
    band_projection_class = factory.get_band_projection_class(
        base_path=base_path)
    assert isinstance(band_projection_class, Procar)


def test_get_potential_class(file_path):
    """
    Test get potential class function
    """
    base_path = file_path("/H/")
    factory = Vasp()
    potential_class = factory.get_potential_class(base_path=base_path)
    assert isinstance(potential_class, Potcar)


def test_get_fermi_energy(file_path):
    """
    Test get fermi energy function
    """
    base_path = file_path("/gec-2d/")
    factory = Vasp()
    fermi_energy = factory.get_fermi_energy(base_path=base_path)
    assert np.isclose(-3.17131785, fermi_energy)


def test_get_number_of_bands(file_path):
    """
    Test get number of bands function
    """
    base_path = file_path("/gec-2d/")
    factory = Vasp()
    bands_num = factory.get_number_of_bands(base_path=base_path)
    assert bands_num == 16


def test_get_number_of_kpoints(file_path):
    """
    Test get number of kpoints function
    """
    base_path = file_path("/gec-2d/")
    factory = Vasp()
    kpoints_num = factory.get_number_of_kpoints(base_path=base_path)
    assert kpoints_num == 12


def test_get_eigenvalues(file_path):
    """
    Test get eigenvalues function
    """
    base_path = file_path("/gec-2d/")
    factory = Vasp()
    eigenvalues = factory.get_eigenvalues(base_path=base_path)
    assert np.isclose(eigenvalues[4][3], -5.397776)


def test_get_runner():
    """
    Test get vasp runner class
    """
    command = ['mpirun', '-np', '4', 'vasp-FEB2016']
    factory = Vasp()
    runner = factory.get_runner(command)
    assert isinstance(runner, VaspRunner)


def test_get_nearest_neighbor_distance(file_path):
    """
    Test get nearest neighbor distance function
    """
    base_path = file_path("/sic-2d/")
    factory = Vasp()
    distance = factory.get_nearest_neighbor_distance(ion_index="1",
                                                     base_path=base_path)
    assert np.isclose(distance, 1.78)


def test_get_number_of_equal_neighbors(file_path):
    """
    Test get number of equal neighbors
    """
    base_path = file_path("/sic-2d/")
    factory = Vasp()
    fake_atoms_map = {"1": "Si", "2": "Si"}
    real_atoms_map = {"1": "Si", "2": "C"}

    assert factory.get_number_of_equal_neighbors(fake_atoms_map,
                                                 symbol="Si",
                                                 base_path=base_path) == 1
    assert factory.get_number_of_equal_neighbors(real_atoms_map,
                                                 symbol="Si",
                                                 base_path=base_path) == 0
