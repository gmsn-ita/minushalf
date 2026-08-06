"""
Test pwscf module
The functions in this file obey the following rules:
test_(what the function is meant to do)_(chemical compound)
"""
import numpy as np
from minushalf.softwares.qe.pwscf import PWSCF


def test_pwscf_aln_2d(file_path):
    """
    Check if the parser for the pwscf.xml file
    is catching the right values for all parsed members
    for AlN (2D).
    """
    filename = file_path("/aln-2d/pwscf.xml")
    pwscf = PWSCF(filename)

    atoms = {'1': 'Al', '2': 'Al', '3': 'N', '4': 'N'}
    fermi_energy = 2.29098605
    nkpts = 40
    nbands = 20
    eigenvalue_firstk_firstband = -12.36076976
    eigenvalue_firstk_lastband = 21.32835345
    nearest_neighbor_distance = 1.890927

    for index, symbol in pwscf.atoms_map.items():
        assert symbol == atoms[index]
    assert np.isclose(pwscf.fermi_energy, fermi_energy)
    assert pwscf.number_of_kpoints == nkpts
    assert pwscf.number_of_bands == nbands
    assert np.isclose(pwscf.eigenvalues[1][0], eigenvalue_firstk_firstband)
    assert np.isclose(pwscf.eigenvalues[1][-1], eigenvalue_firstk_lastband)
    assert len(pwscf.eigenvalues[1]) == nbands
    assert np.isclose(
        pwscf.nearest_neighbor_distance("1"), nearest_neighbor_distance 
    )
    # --- nearest_neighbor_distance: check all ions have entries ---
    for index in pwscf.atoms_map:
        assert len(pwscf.relative_distances[index]) > 0
    assert pwscf.number_of_equal_neighbors(pwscf.atoms_map, "Al") == 1 
    assert pwscf.number_of_equal_neighbors(pwscf.atoms_map, "N") == 1

def test_pwscf_sic_2d(file_path):
    """
    Check if the parser for the pwscf.xml file
    is catching the right values for all parsed members
    for SiC (2D).
    """
    filename = file_path("/sic-2d/pwscf.xml")
    pwscf = PWSCF(filename)

    atoms = {'1': 'Si', '2': 'C'}
    fermi_energy = 2.92451471
    nkpts = 29
    nbands = 8
    eigenvalue_firstk_firstband = -8.29436907
    eigenvalue_firstk_lastband = 6.89856784
    nearest_neighbor_distance = 2.338973

    for index, symbol in pwscf.atoms_map.items():
        assert symbol == atoms[index]
    assert np.isclose(pwscf.fermi_energy, fermi_energy)
    assert pwscf.number_of_kpoints == nkpts
    assert pwscf.number_of_bands == nbands
    assert np.isclose(pwscf.eigenvalues[1][0], eigenvalue_firstk_firstband)
    assert np.isclose(pwscf.eigenvalues[1][-1], eigenvalue_firstk_lastband)
    assert len(pwscf.eigenvalues[1]) == nbands
    assert np.isclose(
        pwscf.nearest_neighbor_distance("1"), nearest_neighbor_distance 
    )
    # --- nearest_neighbor_distance: check all ions have entries ---
    for index in pwscf.atoms_map:
        assert len(pwscf.relative_distances[index]) > 0
    assert pwscf.number_of_equal_neighbors(pwscf.atoms_map, "Si") == 0 
    assert pwscf.number_of_equal_neighbors(pwscf.atoms_map, "C") == 0

