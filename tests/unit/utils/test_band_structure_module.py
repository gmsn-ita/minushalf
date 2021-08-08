"""
Test procar module
The functions in this file obey the following rules:
test_(what the function is meant to do)_(chemical compound)
"""
import numpy as np
from minushalf.utils import BandStructure
from minushalf.softwares.vasp import Procar, Vasprun, Eigenvalues
from minushalf.softwares import Vasp


def test_is_metal_gan_3d(file_path):
    """
    Confirms that ge GaN 3d is not an metal.
    """
    procar_filename = file_path("/gan-3d/PROCAR")
    eigenval_filename = file_path("/gan-3d/EIGENVAL")
    vasprun_filename = file_path("/gan-3d/vasprun.xml")
    procar = Procar(procar_filename)
    vasprun = Vasprun(vasprun_filename)
    eigenval = Eigenvalues(eigenval_filename)

    band_structure = BandStructure(eigenvalues=eigenval.eigenvalues,
                                   fermi_energy=vasprun.fermi_energy,
                                   atoms_map=vasprun.atoms_map,
                                   num_bands=procar.num_bands,
                                   band_projection=procar)

    assert band_structure.is_metal() is False


def test_create_bandstrutcture_gan_3d(file_path):
    """
    Tests create method for GaN 3d
    """
    software_module = Vasp()

    band_structure = BandStructure.create(software_module,
                                          file_path("/gan-3d/"))

    assert isinstance(band_structure, BandStructure)
    assert band_structure.is_metal() is False


def test_band_gap_gan_3d(file_path):
    """
    Tests GaN 3d band gap .
    """
    procar_filename = file_path("/gan-3d/PROCAR")
    eigenval_filename = file_path("/gan-3d/EIGENVAL")
    vasprun_filename = file_path("/gan-3d/vasprun.xml")

    procar = Procar(procar_filename)
    vasprun = Vasprun(vasprun_filename)
    eigenval = Eigenvalues(eigenval_filename)

    band_structure = BandStructure(eigenvalues=eigenval.eigenvalues,
                                   fermi_energy=vasprun.fermi_energy,
                                   atoms_map=vasprun.atoms_map,
                                   num_bands=procar.num_bands,
                                   band_projection=procar)

    assert np.isclose(band_structure.band_gap()["gap"], 1.5380389999999995)



def test_vbm_index_gan_3d(file_path):
    """
    Test if the index of the valence maximum band is correct
    """
    procar_filename = file_path("/gan-3d/PROCAR")
    eigenval_filename = file_path("/gan-3d/EIGENVAL")
    vasprun_filename = file_path("/gan-3d/vasprun.xml")
    kpoint_vbm = 1
    band_vbm = 9

    procar = Procar(procar_filename)
    vasprun = Vasprun(vasprun_filename)
    eigenval = Eigenvalues(eigenval_filename)

    band_structure = BandStructure(eigenvalues=eigenval.eigenvalues,
                                   fermi_energy=vasprun.fermi_energy,
                                   atoms_map=vasprun.atoms_map,
                                   num_bands=procar.num_bands,
                                   band_projection=procar)

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

    procar = Procar(procar_filename)
    vasprun = Vasprun(vasprun_filename)
    eigenval = Eigenvalues(eigenval_filename)

    band_structure = BandStructure(eigenvalues=eigenval.eigenvalues,
                                   fermi_energy=vasprun.fermi_energy,
                                   atoms_map=vasprun.atoms_map,
                                   num_bands=procar.num_bands,
                                   band_projection=procar)

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

    procar = Procar(procar_filename)
    vasprun = Vasprun(vasprun_filename)
    eigenval = Eigenvalues(eigenval_filename)

    band_structure = BandStructure(eigenvalues=eigenval.eigenvalues,
                                   fermi_energy=vasprun.fermi_energy,
                                   atoms_map=vasprun.atoms_map,
                                   num_bands=procar.num_bands,
                                   band_projection=procar)

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

    procar = Procar(procar_filename)
    vasprun = Vasprun(vasprun_filename)
    eigenval = Eigenvalues(eigenval_filename)

    band_structure = BandStructure(eigenvalues=eigenval.eigenvalues,
                                   fermi_energy=vasprun.fermi_energy,
                                   atoms_map=vasprun.atoms_map,
                                   num_bands=procar.num_bands,
                                   band_projection=procar)

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

    procar = Procar(procar_filename)
    vasprun = Vasprun(vasprun_filename)
    eigenval = Eigenvalues(eigenval_filename)

    band_structure = BandStructure(eigenvalues=eigenval.eigenvalues,
                                   fermi_energy=vasprun.fermi_energy,
                                   atoms_map=vasprun.atoms_map,
                                   num_bands=procar.num_bands,
                                   band_projection=procar)

    kpt_9_band_11_projection = {
        "Ga": [0.375, 0.000, 0.000, 0.001, 0.003, 0.007, 0.011, 0.003, 0.034],
        "N": [0.001, 0.080, 0.080, 0.132, 0.000, 0.000, 0.000, 0.000, 0.000]
    }

    band_projection = band_structure.band_projection(9, 11)

    for atom_index, projections in band_projection.items():
        for index, element in enumerate(projections):
            assert np.isclose(element,
                              kpt_9_band_11_projection[atom_index][index])


def test_is_metal_bn_2d(file_path):
    """
    Confirms that  BN 2d is not an metal.
    """
    procar_filename = file_path("/bn-2d/PROCAR")
    eigenval_filename = file_path("/bn-2d/EIGENVAL")
    vasprun_filename = file_path("/bn-2d/vasprun.xml")

    procar = Procar(procar_filename)
    vasprun = Vasprun(vasprun_filename)
    eigenval = Eigenvalues(eigenval_filename)

    band_structure = BandStructure(eigenvalues=eigenval.eigenvalues,
                                   fermi_energy=vasprun.fermi_energy,
                                   atoms_map=vasprun.atoms_map,
                                   num_bands=procar.num_bands,
                                   band_projection=procar)

    assert band_structure.is_metal() is False


def test_band_gap_bn_2d(file_path):
    """
    Test BN 2d band gap.
    """
    procar_filename = file_path("/bn-2d/PROCAR")
    eigenval_filename = file_path("/bn-2d/EIGENVAL")
    vasprun_filename = file_path("/bn-2d/vasprun.xml")

    procar = Procar(procar_filename)
    vasprun = Vasprun(vasprun_filename)
    eigenval = Eigenvalues(eigenval_filename)

    band_structure = BandStructure(eigenvalues=eigenval.eigenvalues,
                                   fermi_energy=vasprun.fermi_energy,
                                   atoms_map=vasprun.atoms_map,
                                   num_bands=procar.num_bands,
                                   band_projection=procar)

    assert np.isclose(band_structure.band_gap()["gap"], 4.756284)


def test_vbm_index_bn_2d(file_path):
    """
    Test if the index of the valence maximum band is correct
    """
    procar_filename = file_path("/bn-2d/PROCAR")
    eigenval_filename = file_path("/bn-2d/EIGENVAL")
    vasprun_filename = file_path("/bn-2d/vasprun.xml")
    kpoint_vbm = 24
    band_vbm = 4

    procar = Procar(procar_filename)
    vasprun = Vasprun(vasprun_filename)
    eigenval = Eigenvalues(eigenval_filename)

    band_structure = BandStructure(eigenvalues=eigenval.eigenvalues,
                                   fermi_energy=vasprun.fermi_energy,
                                   atoms_map=vasprun.atoms_map,
                                   num_bands=procar.num_bands,
                                   band_projection=procar)

    vbm_index = band_structure.vbm_index()
    assert vbm_index[0] == kpoint_vbm
    assert vbm_index[1] == band_vbm

def test_create_bandstrutcture_bn_2d(file_path):
    """
    Tests create method for BN 2d
    """
    softare_module = Vasp()
    kpoint_vbm = 24
    band_vbm = 4

    band_structure = BandStructure.create(softare_module,
                                          file_path("/bn-2d/"))

    vbm_index = band_structure.vbm_index()
    assert isinstance(band_structure, BandStructure)
    assert vbm_index[0] == kpoint_vbm
    assert vbm_index[1] == band_vbm


def test_cbm_index_bn_2d(file_path):
    """
    Test if the index of the conduction minimum band is correct
    """
    procar_filename = file_path("/bn-2d/PROCAR")
    eigenval_filename = file_path("/bn-2d/EIGENVAL")
    vasprun_filename = file_path("/bn-2d/vasprun.xml")
    kpoint_cbm = 1
    band_cbm = 5

    procar = Procar(procar_filename)
    vasprun = Vasprun(vasprun_filename)
    eigenval = Eigenvalues(eigenval_filename)

    band_structure = BandStructure(eigenvalues=eigenval.eigenvalues,
                                   fermi_energy=vasprun.fermi_energy,
                                   atoms_map=vasprun.atoms_map,
                                   num_bands=procar.num_bands,
                                   band_projection=procar)

    cbm_index = band_structure.cbm_index()
    assert cbm_index[0] == kpoint_cbm
    assert cbm_index[1] == band_cbm


def test_vbm_projection_bn_2d(file_path):
    """
    Test if the projection of the valence maximum band is correct
    """
    procar_filename = file_path("/bn-2d/PROCAR")
    eigenval_filename = file_path("/bn-2d/EIGENVAL")
    vasprun_filename = file_path("/bn-2d/vasprun.xml")
    projection = {
        "B": [0.000, 0.000, 0.008, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000],
        "N": [0.000, 0.000, 0.606, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000]
    }

    procar = Procar(procar_filename)
    vasprun = Vasprun(vasprun_filename)
    eigenval = Eigenvalues(eigenval_filename)

    band_structure = BandStructure(eigenvalues=eigenval.eigenvalues,
                                   fermi_energy=vasprun.fermi_energy,
                                   atoms_map=vasprun.atoms_map,
                                   num_bands=procar.num_bands,
                                   band_projection=procar)

    vbm_projection = band_structure.vbm_projection()

    for atom_index, projections in vbm_projection.items():
        for index, element in enumerate(projections):
            assert np.isclose(element, projection[atom_index][index])


def test_cbm_projection_bn_2d(file_path):
    """
    Test if the projection of the valence maximum band is correct
    """
    procar_filename = file_path("/bn-2d/PROCAR")
    eigenval_filename = file_path("/bn-2d/EIGENVAL")
    vasprun_filename = file_path("/bn-2d/vasprun.xml")
    projection = {
        "B": [0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000],
        "N": [0.041, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000]
    }

    procar = Procar(procar_filename)
    vasprun = Vasprun(vasprun_filename)
    eigenval = Eigenvalues(eigenval_filename)

    band_structure = BandStructure(eigenvalues=eigenval.eigenvalues,
                                   fermi_energy=vasprun.fermi_energy,
                                   atoms_map=vasprun.atoms_map,
                                   num_bands=procar.num_bands,
                                   band_projection=procar)

    cbm_projection = band_structure.cbm_projection()

    for atom_index, projections in cbm_projection.items():
        for index, element in enumerate(projections):
            assert np.isclose(element, projection[atom_index][index])


def test_get_band_projection_kpt_9_band_8_bn_2d(file_path):
    """
    Verify if the informations returned about projection
    in the 11ª band of the 9º kpoint is correct for the 2d GaN.
    """
    procar_filename = file_path("/bn-2d/PROCAR")
    eigenval_filename = file_path("/bn-2d/EIGENVAL")
    vasprun_filename = file_path("/bn-2d/vasprun.xml")

    procar = Procar(procar_filename)
    vasprun = Vasprun(vasprun_filename)
    eigenval = Eigenvalues(eigenval_filename)

    band_structure = BandStructure(eigenvalues=eigenval.eigenvalues,
                                   fermi_energy=vasprun.fermi_energy,
                                   atoms_map=vasprun.atoms_map,
                                   num_bands=procar.num_bands,
                                   band_projection=procar)

    kpt_9_band_8_projection = {
        "B": [0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000],
        "N": [0.000, 0.000, 0.032, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000]
    }

    band_projection = band_structure.band_projection(9, 8)

    for atom_index, projections in band_projection.items():
        for index, element in enumerate(projections):
            assert np.isclose(element,
                              kpt_9_band_8_projection[atom_index][index])


def test_is_metal_sic_2d(file_path):
    """
    Confirms that ge SiC 2d is not an metal.
    """
    procar_filename = file_path("/sic-2d/PROCAR")
    eigenval_filename = file_path("/sic-2d/EIGENVAL")
    vasprun_filename = file_path("/sic-2d/vasprun.xml")

    procar = Procar(procar_filename)
    vasprun = Vasprun(vasprun_filename)
    eigenval = Eigenvalues(eigenval_filename)

    band_structure = BandStructure(eigenvalues=eigenval.eigenvalues,
                                   fermi_energy=vasprun.fermi_energy,
                                   atoms_map=vasprun.atoms_map,
                                   num_bands=procar.num_bands,
                                   band_projection=procar)

    assert band_structure.is_metal() is False


def test_band_gap_sic_2d(file_path):
    """
    Test SiC 2d band gap.
    """
    procar_filename = file_path("/sic-2d/PROCAR")
    eigenval_filename = file_path("/sic-2d/EIGENVAL")
    vasprun_filename = file_path("/sic-2d/vasprun.xml")

    procar = Procar(procar_filename)
    vasprun = Vasprun(vasprun_filename)
    eigenval = Eigenvalues(eigenval_filename)

    band_structure = BandStructure(eigenvalues=eigenval.eigenvalues,
                                   fermi_energy=vasprun.fermi_energy,
                                   atoms_map=vasprun.atoms_map,
                                   num_bands=procar.num_bands,
                                   band_projection=procar)

    assert np.isclose(band_structure.band_gap()["gap"], 2.61388)


def test_vbm_index_sic_2d(file_path):
    """
    Test if the index of the valence maximum band is correct
    """
    procar_filename = file_path("/sic-2d/PROCAR")
    eigenval_filename = file_path("/sic-2d/EIGENVAL")
    vasprun_filename = file_path("/sic-2d/vasprun.xml")
    kpoint_vbm = 27
    band_vbm = 4

    procar = Procar(procar_filename)
    vasprun = Vasprun(vasprun_filename)
    eigenval = Eigenvalues(eigenval_filename)

    band_structure = BandStructure(eigenvalues=eigenval.eigenvalues,
                                   fermi_energy=vasprun.fermi_energy,
                                   atoms_map=vasprun.atoms_map,
                                   num_bands=procar.num_bands,
                                   band_projection=procar)

    vbm_index = band_structure.vbm_index()
    assert vbm_index[0] == kpoint_vbm
    assert vbm_index[1] == band_vbm


def test_cbm_index_sic_2d(file_path):
    """
    Test if the index of the conduction minimum band is correct
    """
    procar_filename = file_path("/sic-2d/PROCAR")
    eigenval_filename = file_path("/sic-2d/EIGENVAL")
    vasprun_filename = file_path("/sic-2d/vasprun.xml")
    kpoint_cbm = 15
    band_cbm = 5

    procar = Procar(procar_filename)
    vasprun = Vasprun(vasprun_filename)
    eigenval = Eigenvalues(eigenval_filename)

    band_structure = BandStructure(eigenvalues=eigenval.eigenvalues,
                                   fermi_energy=vasprun.fermi_energy,
                                   atoms_map=vasprun.atoms_map,
                                   num_bands=procar.num_bands,
                                   band_projection=procar)

    cbm_index = band_structure.cbm_index()
    assert cbm_index[0] == kpoint_cbm
    assert cbm_index[1] == band_cbm


def test_vbm_projection_sic_2d(file_path):
    """
    Test if the projection of the valence maximum band is correct
    """
    procar_filename = file_path("/sic-2d/PROCAR")
    eigenval_filename = file_path("/sic-2d/EIGENVAL")
    vasprun_filename = file_path("/sic-2d/vasprun.xml")
    projection = {
        "Si": [0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000],
        "C": [0.000, 0.000, 0.392, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000]
    }

    procar = Procar(procar_filename)
    vasprun = Vasprun(vasprun_filename)
    eigenval = Eigenvalues(eigenval_filename)

    band_structure = BandStructure(eigenvalues=eigenval.eigenvalues,
                                   fermi_energy=vasprun.fermi_energy,
                                   atoms_map=vasprun.atoms_map,
                                   num_bands=procar.num_bands,
                                   band_projection=procar)

    vbm_projection = band_structure.vbm_projection()

    for atom_index, projections in vbm_projection.items():
        for index, element in enumerate(projections):
            assert np.isclose(element, projection[atom_index][index])


def test_cbm_projection_sic_2d(file_path):
    """
    Test if the projection of the valence maximum band is correct
    """
    procar_filename = file_path("/sic-2d/PROCAR")
    eigenval_filename = file_path("/sic-2d/EIGENVAL")
    vasprun_filename = file_path("/sic-2d/vasprun.xml")
    projection = {
        "Si": [0.000, 0.000, 0.254, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000],
        "C": [0.000, 0.000, 0.040, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000]
    }

    procar = Procar(procar_filename)
    vasprun = Vasprun(vasprun_filename)
    eigenval = Eigenvalues(eigenval_filename)

    band_structure = BandStructure(eigenvalues=eigenval.eigenvalues,
                                   fermi_energy=vasprun.fermi_energy,
                                   atoms_map=vasprun.atoms_map,
                                   num_bands=procar.num_bands,
                                   band_projection=procar)

    cbm_projection = band_structure.cbm_projection()

    for atom_index, projections in cbm_projection.items():
        for index, element in enumerate(projections):
            assert np.isclose(element, projection[atom_index][index])


def test_get_band_projection_kpt_9_band_11_sic_2d(file_path):
    """
    Verify if the informations returned about projection
    in the 11ª band of the 9º kpoint is correct for the 2d SiC.
    """
    procar_filename = file_path("/sic-2d/PROCAR")
    eigenval_filename = file_path("/sic-2d/EIGENVAL")
    vasprun_filename = file_path("/sic-2d/vasprun.xml")

    procar = Procar(procar_filename)
    vasprun = Vasprun(vasprun_filename)
    eigenval = Eigenvalues(eigenval_filename)

    band_structure = BandStructure(eigenvalues=eigenval.eigenvalues,
                                   fermi_energy=vasprun.fermi_energy,
                                   atoms_map=vasprun.atoms_map,
                                   num_bands=procar.num_bands,
                                   band_projection=procar)

    kpt_9_band_11_projection = {
        "Si": [0.051, 0.056, 0.000, 0.058, 0.000, 0.000, 0.000, 0.000, 0.000],
        "C": [0.089, 0.002, 0.000, 0.001, 0.000, 0.000, 0.000, 0.000, 0.000]
    }

    band_projection = band_structure.band_projection(9, 11)

    for atom_index, projections in band_projection.items():
        for index, element in enumerate(projections):
            assert np.isclose(element,
                              kpt_9_band_11_projection[atom_index][index])


def test_is_metal_gec_2d(file_path):
    """
    Confirms that  GeC 2d is not an metal.
    """
    procar_filename = file_path("/gec-2d/PROCAR")
    eigenval_filename = file_path("/gec-2d/EIGENVAL")
    vasprun_filename = file_path("/gec-2d/vasprun.xml")

    procar = Procar(procar_filename)
    vasprun = Vasprun(vasprun_filename)
    eigenval = Eigenvalues(eigenval_filename)

    band_structure = BandStructure(eigenvalues=eigenval.eigenvalues,
                                   fermi_energy=vasprun.fermi_energy,
                                   atoms_map=vasprun.atoms_map,
                                   num_bands=procar.num_bands,
                                   band_projection=procar)

    assert band_structure.is_metal() is False


def test_band_gap_gec_2d(file_path):
    """
    Test GeC 2d band gap.
    """
    procar_filename = file_path("/gec-2d/PROCAR")
    eigenval_filename = file_path("/gec-2d/EIGENVAL")
    vasprun_filename = file_path("/gec-2d/vasprun.xml")

    procar = Procar(procar_filename)
    vasprun = Vasprun(vasprun_filename)
    eigenval = Eigenvalues(eigenval_filename)

    band_structure = BandStructure(eigenvalues=eigenval.eigenvalues,
                                   fermi_energy=vasprun.fermi_energy,
                                   atoms_map=vasprun.atoms_map,
                                   num_bands=procar.num_bands,
                                   band_projection=procar)

    assert np.isclose(band_structure.band_gap()["gap"], 2.1609610000000004)


def test_vbm_index_gec_2d(file_path):
    """
    Test if the index of the valence maximum band is correct
    """
    procar_filename = file_path("/gec-2d/PROCAR")
    eigenval_filename = file_path("/gec-2d/EIGENVAL")
    vasprun_filename = file_path("/gec-2d/vasprun.xml")
    kpoint_vbm = 12
    band_vbm = 4

    procar = Procar(procar_filename)
    vasprun = Vasprun(vasprun_filename)
    eigenval = Eigenvalues(eigenval_filename)

    band_structure = BandStructure(eigenvalues=eigenval.eigenvalues,
                                   fermi_energy=vasprun.fermi_energy,
                                   atoms_map=vasprun.atoms_map,
                                   num_bands=procar.num_bands,
                                   band_projection=procar)

    vbm_index = band_structure.vbm_index()
    assert vbm_index[0] == kpoint_vbm
    assert vbm_index[1] == band_vbm


def test_cbm_index_gec_2d(file_path):
    """
    Test if the index of the conduction minimum band is correct
    """
    procar_filename = file_path("/gec-2d/PROCAR")
    eigenval_filename = file_path("/gec-2d/EIGENVAL")
    vasprun_filename = file_path("/gec-2d/vasprun.xml")
    kpoint_cbm = 12
    band_cbm = 5

    procar = Procar(procar_filename)
    vasprun = Vasprun(vasprun_filename)
    eigenval = Eigenvalues(eigenval_filename)

    band_structure = BandStructure(eigenvalues=eigenval.eigenvalues,
                                   fermi_energy=vasprun.fermi_energy,
                                   atoms_map=vasprun.atoms_map,
                                   num_bands=procar.num_bands,
                                   band_projection=procar)

    cbm_index = band_structure.cbm_index()
    assert cbm_index[0] == kpoint_cbm
    assert cbm_index[1] == band_cbm


def test_vbm_projection_gec_2d(file_path):
    """
    Test if the projection of the valence maximum band is correct
    """
    procar_filename = file_path("/gec-2d/PROCAR")
    eigenval_filename = file_path("/gec-2d/EIGENVAL")
    vasprun_filename = file_path("/gec-2d/vasprun.xml")
    projection = {
        "Ge": [0.000, 0.000, 0.000, 0.000, 0.000, 0.028, 0.000, 0.028, 0.000],
        "C": [0.000, 0.000, 0.411, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000]
    }

    procar = Procar(procar_filename)
    vasprun = Vasprun(vasprun_filename)
    eigenval = Eigenvalues(eigenval_filename)

    band_structure = BandStructure(eigenvalues=eigenval.eigenvalues,
                                   fermi_energy=vasprun.fermi_energy,
                                   atoms_map=vasprun.atoms_map,
                                   num_bands=procar.num_bands,
                                   band_projection=procar)

    vbm_projection = band_structure.vbm_projection()

    for atom_index, projections in vbm_projection.items():
        for index, element in enumerate(projections):
            assert np.isclose(element, projection[atom_index][index])


def test_cbm_projection_gec_2d(file_path):
    """
    Test if the projection of the valence maximum band is correct
    """
    procar_filename = file_path("/gec-2d/PROCAR")
    eigenval_filename = file_path("/gec-2d/EIGENVAL")
    vasprun_filename = file_path("/gec-2d/vasprun.xml")
    projection = {
        "Ge": [0.000, 0.000, 0.490, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000],
        "C": [0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000]
    }

    procar = Procar(procar_filename)
    vasprun = Vasprun(vasprun_filename)
    eigenval = Eigenvalues(eigenval_filename)

    band_structure = BandStructure(eigenvalues=eigenval.eigenvalues,
                                   fermi_energy=vasprun.fermi_energy,
                                   atoms_map=vasprun.atoms_map,
                                   num_bands=procar.num_bands,
                                   band_projection=procar)

    cbm_projection = band_structure.cbm_projection()

    for atom_index, projections in cbm_projection.items():
        for index, element in enumerate(projections):
            assert np.isclose(element, projection[atom_index][index])


def test_get_band_projection_kpt_9_band_11_gec_2d(file_path):
    """
    Verify if the informations returned about projection
    in the 11ª band of the 9º kpoint is correct for the 2d GeC.
    """
    procar_filename = file_path("/gec-2d/PROCAR")
    eigenval_filename = file_path("/gec-2d/EIGENVAL")
    vasprun_filename = file_path("/gec-2d/vasprun.xml")

    procar = Procar(procar_filename)
    vasprun = Vasprun(vasprun_filename)
    eigenval = Eigenvalues(eigenval_filename)

    band_structure = BandStructure(eigenvalues=eigenval.eigenvalues,
                                   fermi_energy=vasprun.fermi_energy,
                                   atoms_map=vasprun.atoms_map,
                                   num_bands=procar.num_bands,
                                   band_projection=procar)

    kpt_9_band_11_projection = {
        "Ge": [0.001, 0.001, 0.000, 0.001, 0.002, 0.000, 0.006, 0.000, 0.001],
        "C": [0.001, 0.005, 0.000, 0.016, 0.000, 0.000, 0.000, 0.000, 0.000]
    }

    band_projection = band_structure.band_projection(9, 11)

    for atom_index, projections in band_projection.items():
        for index, element in enumerate(projections):
            assert np.isclose(element,
                              kpt_9_band_11_projection[atom_index][index])


def test_is_metal_aln_2d(file_path):
    """
    Confirms that AlN 2d is not an metal.
    """
    procar_filename = file_path("/aln-2d/PROCAR")
    eigenval_filename = file_path("/aln-2d/EIGENVAL")
    vasprun_filename = file_path("/aln-2d/vasprun.xml")

    procar = Procar(procar_filename)
    vasprun = Vasprun(vasprun_filename)
    eigenval = Eigenvalues(eigenval_filename)

    band_structure = BandStructure(eigenvalues=eigenval.eigenvalues,
                                   fermi_energy=vasprun.fermi_energy,
                                   atoms_map=vasprun.atoms_map,
                                   num_bands=procar.num_bands,
                                   band_projection=procar)

    assert band_structure.is_metal() is False


def test_band_gap_aln_2d(file_path):
    """
    Test AlN 2d band gap
    """
    procar_filename = file_path("/aln-2d/PROCAR")
    eigenval_filename = file_path("/aln-2d/EIGENVAL")
    vasprun_filename = file_path("/aln-2d/vasprun.xml")

    procar = Procar(procar_filename)
    vasprun = Vasprun(vasprun_filename)
    eigenval = Eigenvalues(eigenval_filename)

    band_structure = BandStructure(eigenvalues=eigenval.eigenvalues,
                                   fermi_energy=vasprun.fermi_energy,
                                   atoms_map=vasprun.atoms_map,
                                   num_bands=procar.num_bands,
                                   band_projection=procar)

    assert np.isclose(band_structure.band_gap()["gap"], 2.924163)


def test_vbm_index_aln_2d(file_path):
    """
    Test if the index of the valence maximum band is correct
    """
    procar_filename = file_path("/aln-2d/PROCAR")
    eigenval_filename = file_path("/aln-2d/EIGENVAL")
    vasprun_filename = file_path("/aln-2d/vasprun.xml")
    kpoint_vbm = 16
    band_vbm = 4

    procar = Procar(procar_filename)
    vasprun = Vasprun(vasprun_filename)
    eigenval = Eigenvalues(eigenval_filename)

    band_structure = BandStructure(eigenvalues=eigenval.eigenvalues,
                                   fermi_energy=vasprun.fermi_energy,
                                   atoms_map=vasprun.atoms_map,
                                   num_bands=procar.num_bands,
                                   band_projection=procar)

    vbm_index = band_structure.vbm_index()
    assert vbm_index[0] == kpoint_vbm
    assert vbm_index[1] == band_vbm


def test_cbm_index_aln_2d(file_path):
    """
    Test if the index of the conduction minimum band is correct
    """
    procar_filename = file_path("/aln-2d/PROCAR")
    eigenval_filename = file_path("/aln-2d/EIGENVAL")
    vasprun_filename = file_path("/aln-2d/vasprun.xml")
    kpoint_cbm = 1
    band_cbm = 5

    procar = Procar(procar_filename)
    vasprun = Vasprun(vasprun_filename)
    eigenval = Eigenvalues(eigenval_filename)

    band_structure = BandStructure(eigenvalues=eigenval.eigenvalues,
                                   fermi_energy=vasprun.fermi_energy,
                                   atoms_map=vasprun.atoms_map,
                                   num_bands=procar.num_bands,
                                   band_projection=procar)

    cbm_index = band_structure.cbm_index()
    assert cbm_index[0] == kpoint_cbm
    assert cbm_index[1] == band_cbm


def test_vbm_projection_aln_2d(file_path):
    """
    Test if the projection of the valence maximum band is correct
    """
    procar_filename = file_path("/aln-2d/PROCAR")
    eigenval_filename = file_path("/aln-2d/EIGENVAL")
    vasprun_filename = file_path("/aln-2d/vasprun.xml")
    projection = {
        "Al": [0.000, 0.000, 0.002, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000],
        "N": [0.000, 0.000, 0.510, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000]
    }

    procar = Procar(procar_filename)
    vasprun = Vasprun(vasprun_filename)
    eigenval = Eigenvalues(eigenval_filename)

    band_structure = BandStructure(eigenvalues=eigenval.eigenvalues,
                                   fermi_energy=vasprun.fermi_energy,
                                   atoms_map=vasprun.atoms_map,
                                   num_bands=procar.num_bands,
                                   band_projection=procar)

    vbm_projection = band_structure.vbm_projection()

    for atom_index, projections in vbm_projection.items():
        for index, element in enumerate(projections):
            assert np.isclose(element, projection[atom_index][index])


def test_cbm_projection_aln_2d(file_path):
    """
    Test if the projection of the valence maximum band is correct
    """
    procar_filename = file_path("/aln-2d/PROCAR")
    eigenval_filename = file_path("/aln-2d/EIGENVAL")
    vasprun_filename = file_path("/aln-2d/vasprun.xml")
    projection = {
        "Al": [0.055, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000],
        "N": [0.132, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000]
    }

    procar = Procar(procar_filename)
    vasprun = Vasprun(vasprun_filename)
    eigenval = Eigenvalues(eigenval_filename)

    band_structure = BandStructure(eigenvalues=eigenval.eigenvalues,
                                   fermi_energy=vasprun.fermi_energy,
                                   atoms_map=vasprun.atoms_map,
                                   num_bands=procar.num_bands,
                                   band_projection=procar)

    cbm_projection = band_structure.cbm_projection()

    for atom_index, projections in cbm_projection.items():
        for index, element in enumerate(projections):
            assert np.isclose(element, projection[atom_index][index])


def test_get_band_projection_kpt_9_band_8_aln_2d(file_path):
    """
    Verify if the informations returned about projection
    in the 11ª band of the 9º kpoint is correct for the 2d GeC.
    """
    procar_filename = file_path("/aln-2d/PROCAR")
    eigenval_filename = file_path("/aln-2d/EIGENVAL")
    vasprun_filename = file_path("/aln-2d/vasprun.xml")

    procar = Procar(procar_filename)
    vasprun = Vasprun(vasprun_filename)
    eigenval = Eigenvalues(eigenval_filename)

    band_structure = BandStructure(eigenvalues=eigenval.eigenvalues,
                                   fermi_energy=vasprun.fermi_energy,
                                   atoms_map=vasprun.atoms_map,
                                   num_bands=procar.num_bands,
                                   band_projection=procar)

    kpt_9_band_8_projection = {
        "Al": [0.000, 0.000, 0.033, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000],
        "N": [0.000, 0.000, 0.011, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000]
    }

    band_projection = band_structure.band_projection(9, 8)

    for atom_index, projections in band_projection.items():
        for index, element in enumerate(projections):
            assert np.isclose(element,
                              kpt_9_band_8_projection[atom_index][index])
