"""
Test Quantum ESPRESSO factory module
"""
import numpy as np
import os
from minushalf.softwares.qe.potential import Potential
from minushalf.softwares.qe.projoutput import ProjOutput
from minushalf.softwares.qe.pwscf import PWSCF
from minushalf.softwares.qe.runner import QERunner
from minushalf.softwares.qe.qe_factory import QE


def test_get_atoms_map(file_path):
    """
    Test get atoms map function
    """
    base_path = file_path("/aln-2d/")
    factory = QE()
    atoms_map = factory.get_atoms_map(base_path=base_path)
    assert atoms_map["1"] == "Al"
    assert atoms_map["3"] == "N"


def test_get_band_projection_class(file_path):
    """
    Test get band projection class function
    """
    filename = file_path("aln-2d/proj.projwfc_up")
    base_path = file_path("/aln-2d/")
    factory = QE()
    band_projection_class = factory.get_band_projection_class(filename=filename,
        base_path=base_path)
    assert isinstance(band_projection_class, ProjOutput)


def test_get_potential_class(file_path):
    """
    Test get potential class function
    """
    filename = file_path("Al/Al.upf")
    base_path = file_path("/Al/")
    factory = QE()
    potential_class = factory.get_potential_class(filename=filename,base_path=base_path)
    assert isinstance(potential_class, Potential)


def test_get_fermi_energy(file_path):
    """
    Test get fermi energy function
    """
    filename = file_path("/aln-2d/pwscf.xml")
    base_path = file_path("/aln-2d/")
    factory = QE()
    fermi_energy = factory.get_fermi_energy(filename=filename,base_path=base_path)
    assert np.isclose(2.29098605, fermi_energy)


def test_get_number_of_bands(file_path):
    """
    Test get number of bands function
    """
    filename = file_path("/aln-2d/pwscf.xml")
    base_path = file_path("/aln-2d/")
    factory = QE()
    bands_num = factory.get_number_of_bands(filename=filename,base_path=base_path)
    assert bands_num == 20


def test_get_number_of_kpoints(file_path):
    """
    Test get number of kpoints function
    """
    filename = file_path("/aln-2d/pwscf.xml")
    base_path = file_path("/aln-2d/")
    factory = QE()
    kpoints_num = factory.get_number_of_kpoints(filename=filename,base_path=base_path)
    assert kpoints_num == 40


def test_get_eigenvalues(file_path):
    """
    Test get eigenvalues function
    """
    filename = file_path("/aln-2d/pwscf.xml")
    base_path = file_path("/aln-2d/")
    factory = QE()
    eigenvalues = factory.get_eigenvalues(filename=filename,base_path=base_path)
    assert np.isclose(eigenvalues[4][3], -0.526774853)


def test_get_runner():
    """
    Test get vasp runner class
    """
    command = ['mpirun', '-np', '4', 'pw.x']
    factory = QE()
    runner = factory.get_runner(command)
    assert isinstance(runner, QERunner)


def test_get_nearest_neighbor_distance(file_path):
    """
    Test get nearest neighbor distance function
    """
    filename = file_path("/sic-2d/pwscf.xml")
    base_path = file_path("/sic-2d/")
    factory = QE()
    distance = factory.get_nearest_neighbor_distance(ion_index="1",filename=filename,
                                                     base_path=base_path)
    assert np.isclose(distance, 2.338973)


def test_get_number_of_equal_neighbors(file_path):
    """
    Test get number of equal neighbors
    """
    filename = file_path("/sic-2d/pwscf.xml")
    base_path = file_path("/sic-2d/")
    factory = QE()
    fake_atoms_map = {"1": "Si", "2": "Si"}
    real_atoms_map = {"1": "Si", "2": "C"}

    assert factory.get_number_of_equal_neighbors(fake_atoms_map,
                                                 symbol="Si",
                                                 filename=filename,
                                                 base_path=base_path) == 1
    assert factory.get_number_of_equal_neighbors(real_atoms_map,
                                                 symbol="Si",
                                                 filename=filename,
                                                 base_path=base_path) == 0
