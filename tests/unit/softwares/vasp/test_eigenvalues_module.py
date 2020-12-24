"""
Test Eigenval module
The functions in this file obey the following rules:
test_(what the function is meant to do)_(chemical compound)
"""
from minushalf.softwares.vasp import Eigenvalues


def test_parse_eigenvalues_gan_3d(file_path):
    """
    Test if the eigenvalues module
    collect all the eigenvalues.
    """
    filename = file_path("/gan-3d/EIGENVAL")
    eigenval = Eigenvalues(filename)

    number_of_bands = 16
    number_of_kpoints = 10
    for _, band in eigenval.eigenvalues.items():

        assert len(band) == number_of_bands

    assert len(eigenval.eigenvalues) == number_of_kpoints
