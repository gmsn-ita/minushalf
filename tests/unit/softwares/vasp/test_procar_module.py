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


def test_parse_procar_header_bn_2d(file_path):
    """
    Test if the procar header is correctly parsed.
    The header of poscar has 3 lines and contain informations
    about the total number of kpoints and the number of
    bands for each kpoint.
    """
    filename = file_path("/bn-2d/PROCAR")
    procar = Procar(filename)
    number_of_kpoints = 24
    number_of_bands = 8
    size_procar_header = 3

    assert procar.num_kpoints == number_of_kpoints
    assert procar.num_bands == number_of_bands
    assert procar.size_procar_header == size_procar_header


def test_get_band_projection_kpt_1_band_5_bn_2d(file_path):
    """
    Verify if the informations returned about projection
    in the 5ª band of the 1º kpoint is correct for the 2d BN.
    """
    filename = file_path("/bn-2d/PROCAR")
    procar = Procar(filename)

    kpt_1_band_5_projection = {
        '1': [0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000],
        '2': [0.041, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000]
    }

    band_projection = procar.get_band_projection(1, 5)

    for atom_index, projections in band_projection.items():
        for index, element in enumerate(projections):
            assert np.isclose(element,
                              kpt_1_band_5_projection[atom_index][index])


def test_get_band_projection_kpt_9_band_8_bn_2d(file_path):
    """
    Verify if the informations returned about projection
    in the 11ª band of the 9º kpoint is correct for the 2d GaN.
    """
    filename = file_path("/bn-2d/PROCAR")
    procar = Procar(filename)

    kpt_9_band_8_projection = {
        '1': [0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000],
        '2': [0.000, 0.000, 0.032, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000]
    }

    band_projection = procar.get_band_projection(9, 8)

    for atom_index, projections in band_projection.items():
        for index, element in enumerate(projections):
            assert np.isclose(element,
                              kpt_9_band_8_projection[atom_index][index])


def test_parse_procar_header_sic_2d(file_path):
    """
    Test if the procar header is correctly parsed.
    The header of poscar has 3 lines and contain informations
    about the total number of kpoints and the number of
    bands for each kpoint.
    """
    filename = file_path("/sic-2d/PROCAR")
    procar = Procar(filename)
    number_of_kpoints = 27
    number_of_bands = 16
    size_procar_header = 3

    assert procar.num_kpoints == number_of_kpoints
    assert procar.num_bands == number_of_bands
    assert procar.size_procar_header == size_procar_header


def test_get_band_projection_kpt_1_band_5_sic_2d(file_path):
    """
    Verify if the informations returned about projection
    in the 5ª band of the 1º kpoint is correct for the 2d SiC.
    """
    filename = file_path("/sic-2d/PROCAR")
    procar = Procar(filename)

    kpt_1_band_5_projection = {
        '1': [0.015, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000],
        '2': [0.112, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000]
    }

    band_projection = procar.get_band_projection(1, 5)

    for atom_index, projections in band_projection.items():
        for index, element in enumerate(projections):
            assert np.isclose(element,
                              kpt_1_band_5_projection[atom_index][index])


def test_get_band_projection_kpt_9_band_11_sic_2d(file_path):
    """
    Verify if the informations returned about projection
    in the 11ª band of the 9º kpoint is correct for the 2d SiC.
    """
    filename = file_path("/sic-2d/PROCAR")
    procar = Procar(filename)

    kpt_9_band_11_projection = {
        '1': [0.051, 0.056, 0.000, 0.058, 0.000, 0.000, 0.000, 0.000, 0.000],
        '2': [0.089, 0.002, 0.000, 0.001, 0.000, 0.000, 0.000, 0.000, 0.000]
    }

    band_projection = procar.get_band_projection(9, 11)

    for atom_index, projections in band_projection.items():
        for index, element in enumerate(projections):
            assert np.isclose(element,
                              kpt_9_band_11_projection[atom_index][index])


def test_parse_procar_header_gec_2d(file_path):
    """
    Test if the procar header is correctly parsed.
    The header of poscar has 3 lines and contain informations
    about the total number of kpoints and the number of
    bands for each kpoint.
    """
    filename = file_path("/gec-2d/PROCAR")
    procar = Procar(filename)
    number_of_kpoints = 12
    number_of_bands = 16
    size_procar_header = 3

    assert procar.num_kpoints == number_of_kpoints
    assert procar.num_bands == number_of_bands
    assert procar.size_procar_header == size_procar_header


def test_get_band_projection_kpt_1_band_5_gec_2d(file_path):
    """
    Verify if the informations returned about projection
    in the 5ª band of the 1º kpoint is correct for the 2d GeC.
    """
    filename = file_path("/gec-2d/PROCAR")
    procar = Procar(filename)

    kpt_1_band_5_projection = {
        '1': [0.098, 0.000, 0.000, 0.000, 0.000, 0.000, 0.060, 0.000, 0.000],
        '2': [0.244, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000]
    }

    band_projection = procar.get_band_projection(1, 5)

    for atom_index, projections in band_projection.items():
        for index, element in enumerate(projections):
            assert np.isclose(element,
                              kpt_1_band_5_projection[atom_index][index])


def test_get_band_projection_kpt_9_band_11_gec_2d(file_path):
    """
    Verify if the informations returned about projection
    in the 11ª band of the 9º kpoint is correct for the 2d SiC.
    """
    filename = file_path("/gec-2d/PROCAR")
    procar = Procar(filename)

    kpt_9_band_11_projection = {
        '1': [0.001, 0.001, 0.000, 0.001, 0.002, 0.000, 0.006, 0.000, 0.001],
        '2': [0.001, 0.005, 0.000, 0.016, 0.000, 0.000, 0.000, 0.000, 0.000]
    }

    band_projection = procar.get_band_projection(9, 11)

    for atom_index, projections in band_projection.items():
        for index, element in enumerate(projections):
            assert np.isclose(element,
                              kpt_9_band_11_projection[atom_index][index])


def test_parse_procar_header_aln_2d(file_path):
    """
    Test if the procar header is correctly parsed.
    The header of poscar has 3 lines and contain informations
    about the total number of kpoints and the number of
    bands for each kpoint.
    """
    filename = file_path("/aln-2d/PROCAR")
    procar = Procar(filename)
    number_of_kpoints = 16
    number_of_bands = 8
    size_procar_header = 3

    assert procar.num_kpoints == number_of_kpoints
    assert procar.num_bands == number_of_bands
    assert procar.size_procar_header == size_procar_header


def test_get_band_projection_kpt_1_band_5_aln_2d(file_path):
    """
    Verify if the informations returned about projection
    in the 5ª band of the 1º kpoint is correct for the 2d AlN.
    """
    filename = file_path("/aln-2d/PROCAR")
    procar = Procar(filename)

    kpt_1_band_5_projection = {
        '1': [0.055, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000],
        '2': [0.132, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000]
    }

    band_projection = procar.get_band_projection(1, 5)

    for atom_index, projections in band_projection.items():
        for index, element in enumerate(projections):
            assert np.isclose(element,
                              kpt_1_band_5_projection[atom_index][index])


def test_get_band_projection_kpt_9_band_8_bn_2d(file_path):
    """
    Verify if the informations returned about projection
    in the 11ª band of the 9º kpoint is correct for the 2d GaN.
    """
    filename = file_path("/aln-2d/PROCAR")
    procar = Procar(filename)

    kpt_9_band_8_projection = {
        '1': [0.000, 0.000, 0.033, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000],
        '2': [0.000, 0.000, 0.011, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000]
    }

    band_projection = procar.get_band_projection(9, 8)

    for atom_index, projections in band_projection.items():
        for index, element in enumerate(projections):
            assert np.isclose(element,
                              kpt_9_band_8_projection[atom_index][index])