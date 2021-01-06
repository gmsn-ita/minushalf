"""
Test vasp factory module
"""
import numpy as np
from minushalf.softwares.vasp import Procar, Potcar, VaspRunner
from minushalf.softwares import VaspFactory


def test_get_atoms_map(file_path):
    """
    Test get atoms map function
    """
    base_path = file_path("/gec-2d/")
    factory = VaspFactory()
    atoms_map = factory.get_atoms_map(base_path=base_path)
    assert atoms_map["1"] == "Ge"
    assert atoms_map["2"] == "C"


def test_get_band_projection_class(file_path):
    """
    Test get band projection class function
    """
    base_path = file_path("/gec-2d/")
    factory = VaspFactory()
    band_projection_class = factory.get_band_projection_class(
        base_path=base_path)
    assert isinstance(band_projection_class, Procar)


def test_get_potential_class(file_path):
    """
    Test get potential class function
    """
    base_path = file_path("/H/")
    factory = VaspFactory()
    potential_class = factory.get_potential_class(base_path=base_path)
    assert isinstance(potential_class, Potcar)


def test_get_fermi_energy(file_path):
    """
    Test get fermi energy function
    """
    base_path = file_path("/gec-2d/")
    factory = VaspFactory()
    fermi_energy = factory.get_fermi_energy(base_path=base_path)
    assert np.isclose(-3.17131785, fermi_energy)


def test_get_number_of_bands(file_path):
    """
    Test get number of bands function
    """
    base_path = file_path("/gec-2d/")
    factory = VaspFactory()
    bands_num = factory.get_number_of_bands(base_path=base_path)
    assert bands_num == 16


def test_get_number_of_kpoints(file_path):
    """
    Test get number of kpoints function
    """
    base_path = file_path("/gec-2d/")
    factory = VaspFactory()
    kpoints_num = factory.get_number_of_kpoints(base_path=base_path)
    assert kpoints_num == 12


def test_get_eigenvalues(file_path):
    """
    Test get eigenvalues function
    """
    base_path = file_path("/gec-2d/")
    factory = VaspFactory()
    eigenvalues = factory.get_eigenvalues(base_path=base_path)
    assert np.isclose(eigenvalues[4][3], -5.397776)


def test_get_runner():
    """
    Test get vasp rnner class
    """
    factory = VaspFactory()
    runner = factory.get_runner()
    assert isinstance(runner, VaspRunner)
