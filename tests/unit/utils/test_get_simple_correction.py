"""
Test get simple correction function
"""
from minushalf.softwares.vasp import (Procar, Vasprun, Eigenvalues)
from minushalf.utils import (get_fractionary_correction_indexes,
                             projection_to_df, BandStructure)


def test_aln_2d_vbm(file_path):
    """
    Test with AlN-2d
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

    vbm_projection = band_structure.vbm_projection()
    vbm_df = projection_to_df(vbm_projection)
    correction_indexes = get_fractionary_correction_indexes(vbm_df)
    assert correction_indexes["N"][0] == "p"


def test_aln_2d_cbm(file_path):
    """
    Test with AlN-2d
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

    cbm_projection = band_structure.cbm_projection()
    cbm_df = projection_to_df(cbm_projection)
    correction_indexes = get_fractionary_correction_indexes(cbm_df)
    assert correction_indexes["N"][0] == "s"


def test_gec_2d_vbm(file_path):
    """
    Test with GeC-2d
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

    vbm_projection = band_structure.vbm_projection()
    vbm_df = projection_to_df(vbm_projection)
    correction_indexes = get_fractionary_correction_indexes(vbm_df)
    assert correction_indexes["C"][0] == "p"


def test_gec_2d_cbm(file_path):
    """
    Test with GeC-2d
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

    cbm_projection = band_structure.cbm_projection()
    cbm_df = projection_to_df(cbm_projection)
    correction_indexes = get_fractionary_correction_indexes(cbm_df)
    assert correction_indexes["Ge"][0] == "p"


def test_bn_2d_vbm(file_path):
    """
    Test with BN-2d
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

    vbm_projection = band_structure.vbm_projection()
    vbm_df = projection_to_df(vbm_projection)
    correction_indexes = get_fractionary_correction_indexes(vbm_df)
    assert correction_indexes["N"][0] == "p"


def test_bn_2d_cbm(file_path):
    """
    Test with BN-2d
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

    cbm_projection = band_structure.cbm_projection()
    cbm_df = projection_to_df(cbm_projection)
    correction_indexes = get_fractionary_correction_indexes(cbm_df)
    assert correction_indexes["N"][0] == "s"
