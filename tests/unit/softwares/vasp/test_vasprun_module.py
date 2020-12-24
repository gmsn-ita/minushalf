"""
Test vasprun module
The functions in this file obey the following rules:
test_(what the function is meant to do)_(chemical compound)
"""
import numpy as np
from minushalf.softwares.vasp import Vasprun


def test_vasprun_parser_for_fermi_energy_and_atoms_gan_3d(file_path):
    """
    Check if the parser for the vasprun.xml file
    is catching the right values for atom informations
    and fermi energy.
    """
    filename = file_path("/gan-3d/vasprun.xml")
    vasprun = Vasprun(filename)

    fermi_energy = 5.06822674
    atoms = {"1": "Ga", "2": "N"}

    assert np.isclose(fermi_energy, vasprun.fermi_energy)
    for index, symbol in vasprun.atoms_map.items():
        assert symbol == atoms[index]


def test_vasprun_parser_for_fermi_energy_and_atoms_gan_2d(file_path):
    """
    Check if the parser for the vasprun.xml file
    is catching the right values for atom informations
    and fermi energy.
    """
    filename = file_path("/gan-2d/vasprun.xml")
    vasprun = Vasprun(filename)

    fermi_energy = -6.12961499
    atoms = {"1": "Ga", "2": "N"}

    assert np.isclose(fermi_energy, vasprun.fermi_energy)
    for index, symbol in vasprun.atoms_map.items():
        assert symbol == atoms[index]


def test_vasprun_parser_for_fermi_energy_and_atoms_sic_2d(file_path):
    """
    Check if the parser for the vasprun.xml file
    is catching the right values for atom informations
    and fermi energy.
    """
    filename = file_path("/sic-2d/vasprun.xml")
    vasprun = Vasprun(filename)

    fermi_energy = -3.10704937
    atoms = {"1": "Si", "2": "C"}

    assert np.isclose(fermi_energy, vasprun.fermi_energy)
    for index, symbol in vasprun.atoms_map.items():
        assert symbol == atoms[index]


def test_vasprun_parser_for_fermi_energy_and_atoms_gec_2d(file_path):
    """
    Check if the parser for the vasprun.xml file
    is catching the right values for atom informations
    and fermi energy.
    """
    filename = file_path("/gec-2d/vasprun.xml")
    vasprun = Vasprun(filename)

    fermi_energy = -3.17131785
    atoms = {"1": "Ge", "2": "C"}

    assert np.isclose(fermi_energy, vasprun.fermi_energy)
    for index, symbol in vasprun.atoms_map.items():
        assert symbol == atoms[index]
