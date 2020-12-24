"""
Test procar module
The functions in this file obey the following rules:
test_(what the function is meant to do)_(chemical compound)
"""
import numpy as np
from minushalf.softwares.vasp import Procar


def test_parse_procar_header_gan_3d(file_path):
    """
    Test if the procar header is correctly parsed.
    The header of poscar has 3 lines and contain informations
    about the total number of kpoints and the number of
    bands for each kpoint.
    """
    filename = file_path("/gan-3d/PROCAR")
    procar = Procar(filename)
    number_of_kpoints = 10
    number_of_bands = 16
    size_procar_header = 3

    assert procar.num_kpoints == number_of_kpoints
    assert procar.num_bands == number_of_bands
    assert procar.size_procar_header == size_procar_header


def test_get_band_projection_kpt_1_band_5_gan_3d(file_path):
    """
    Verify if the informations returned about projection
    in the 5ª band of the 1º kpoint is correct for the 3d GaN.
    """
    filename = file_path("/gan-3d/PROCAR")
    procar = Procar(filename)

    kpt_1_band_5_projection = {
        '1': [0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.249, 0.000, 0.743],
        '2': [0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000]
    }

    band_projection = procar.get_band_projection(1, 5)

    for atom_index, projections in band_projection.items():
        for index, element in enumerate(projections):
            assert np.isclose(element,
                              kpt_1_band_5_projection[atom_index][index])


def test_get_band_projection_kpt_9_band_11_gan_3d(file_path):
    """
    Verify if the informations returned about projection
    in the 11ª band of the 9º kpoint is correct for the 3d GaN.
    """
    filename = file_path("/gan-3d/PROCAR")
    procar = Procar(filename)

    kpt_9_band_11_projection = {
        '1': [0.375, 0.000, 0.000, 0.001, 0.003, 0.007, 0.011, 0.003, 0.034],
        '2': [0.001, 0.080, 0.080, 0.132, 0.000, 0.000, 0.000, 0.000, 0.000]
    }

    band_projection = procar.get_band_projection(9, 11)

    for atom_index, projections in band_projection.items():
        for index, element in enumerate(projections):
            assert np.isclose(element,
                              kpt_9_band_11_projection[atom_index][index])
