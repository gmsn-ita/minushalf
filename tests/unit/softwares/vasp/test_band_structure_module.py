"""
Test procar module
The functions in this file obey the following rules:
test_(what the function is meant to do)_(chemical compound)
"""
import numpy as np
from minushalf.softwares.vasp import Procar, Vasprun, Eigenvalues, BandStructure


def test_is_metal_gan_3d(file_path):
    """
    Confirms that ge GaN 3d is not an metal.
    """
    procar_filename = file_path("/gan-3d/PROCAR")
    eigenval_filename = file_path("/gan-3d/EIGENVAL")
    vasprun_filename = file_path("/gan-3d/vasprun.xml")

    band_structure = BandStructure(procar=Procar(procar_filename),
                                   vasprun=Vasprun(vasprun_filename),
                                   eigenval=Eigenvalues(eigenval_filename))

    assert band_structure.is_metal() == False


def test_vbm_index_gan_3d(file_path):
    """
    Test if the index of the valence maximum band is correct
    """
    procar_filename = file_path("/gan-3d/PROCAR")
    eigenval_filename = file_path("/gan-3d/EIGENVAL")
    vasprun_filename = file_path("/gan-3d/vasprun.xml")
    kpoint_vbm = 1
    band_vbm = 9

    band_structure = BandStructure(procar=Procar(procar_filename),
                                   vasprun=Vasprun(vasprun_filename),
                                   eigenval=Eigenvalues(eigenval_filename))

    vbm_index = band_structure.vbm_index()
    assert vbm_index[0] == kpoint_vbm
    assert vbm_index[1] == band_vbm


def test_cbm_index_gan_3d(file_path):
    """
    Test if the index of the conduction minimum band is correct
    """
    procar_filename = file_path("/gan-3d/PROCAR")
    eigenval_filename = file_path("/gan-3d/EIGENVAL")
    vasprun_filename = file_path("/gan-3d/vasprun.xml")
    kpoint_cbm = 1
    band_cbm = 10

    band_structure = BandStructure(procar=Procar(procar_filename),
                                   vasprun=Vasprun(vasprun_filename),
                                   eigenval=Eigenvalues(eigenval_filename))

    cbm_index = band_structure.cbm_index()
    assert cbm_index[0] == kpoint_cbm
    assert cbm_index[1] == band_cbm


def test_vbm_projection_gan_3d(file_path):
    """
    Test if the projection of the valence maximum band is correct
    """
    procar_filename = file_path("/gan-3d/PROCAR")
    eigenval_filename = file_path("/gan-3d/EIGENVAL")
    vasprun_filename = file_path("/gan-3d/vasprun.xml")
    projection = {
        "Ga": [0.000, 0.006, 0.034, 0.006, 0.094, 0.016, 0.000, 0.016, 0.000],
        "N": [0.000, 0.076, 0.446, 0.074, 0.000, 0.000, 0.000, 0.000, 0.000]
    }

    band_structure = BandStructure(procar=Procar(procar_filename),
                                   vasprun=Vasprun(vasprun_filename),
                                   eigenval=Eigenvalues(eigenval_filename))

    vbm_projection = band_structure.vbm_projection()

    for atom_index, projections in vbm_projection.items():
        for index, element in enumerate(projections):
            assert np.isclose(element, projection[atom_index][index])


def test_cbm_projection_gan_3d(file_path):
    """
    Test if the projection of the valence maximum band is correct
    """
    procar_filename = file_path("/gan-3d/PROCAR")
    eigenval_filename = file_path("/gan-3d/EIGENVAL")
    vasprun_filename = file_path("/gan-3d/vasprun.xml")
    projection = {
        "Ga": [0.276, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000],
        "N": [0.355, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000]
    }

    band_structure = BandStructure(procar=Procar(procar_filename),
                                   vasprun=Vasprun(vasprun_filename),
                                   eigenval=Eigenvalues(eigenval_filename))

    cbm_projection = band_structure.cbm_projection()

    for atom_index, projections in cbm_projection.items():
        for index, element in enumerate(projections):
            assert np.isclose(element, projection[atom_index][index])


def test_get_band_projection_kpt_9_band_11_gan_3d(file_path):
    """
    Verify if the informations returned about projection
    in the 11ª band of the 9º kpoint is correct for the 3d GaN.
    """
    procar_filename = file_path("/gan-3d/PROCAR")
    eigenval_filename = file_path("/gan-3d/EIGENVAL")
    vasprun_filename = file_path("/gan-3d/vasprun.xml")

    band_structure = BandStructure(procar=Procar(procar_filename),
                                   vasprun=Vasprun(vasprun_filename),
                                   eigenval=Eigenvalues(eigenval_filename))

    kpt_9_band_11_projection = {
        "Ga": [0.375, 0.000, 0.000, 0.001, 0.003, 0.007, 0.011, 0.003, 0.034],
        "N": [0.001, 0.080, 0.080, 0.132, 0.000, 0.000, 0.000, 0.000, 0.000]
    }

    band_projection = band_structure.band_projection(9, 11)

    for atom_index, projections in band_projection.items():
        for index, element in enumerate(projections):
            assert np.isclose(element,
                              kpt_9_band_11_projection[atom_index][index])
