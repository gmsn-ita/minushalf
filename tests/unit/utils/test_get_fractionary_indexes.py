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
    assert correction_indexes["Al"][0] == "s"


def test_aln_2d_cbm_treshold_29(file_path):
    """
    Test with AlN-2d with treshold 29
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
    correction_indexes = get_fractionary_correction_indexes(cbm_df,
                                                            treshold=28)
    assert correction_indexes["N"][0] == "s"
    assert correction_indexes["Al"][0] == "s"


def test_aln_2d_cbm_treshold_30(file_path):
    """
    Test with AlN-2d with treshold 29
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
    correction_indexes = get_fractionary_correction_indexes(cbm_df,
                                                            treshold=29)
    assert correction_indexes["N"][0] == "s"
    assert len(correction_indexes["Al"]) == 0


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
    assert correction_indexes["Ge"][0] == "d"


def test_gec_2d_vbm_changing_treshold_12(file_path):
    """
    Test with GeC-2d with treshold_12
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
    correction_indexes = get_fractionary_correction_indexes(vbm_df,
                                                            treshold=11)
    assert correction_indexes["C"][0] == "p"
    assert correction_indexes["Ge"][0] == "d"


def test_gec_2d_vbm_changing_treshold_13(file_path):
    """
    Test with GeC-2d with treshold_13
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
    correction_indexes = get_fractionary_correction_indexes(vbm_df,
                                                            treshold=12)
    assert correction_indexes["C"][0] == "p"
    assert len(correction_indexes["Ge"]) == 0


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
    assert correction_indexes.get("B") == None


def test_bn_2d_vbm_without_treshold(file_path):
    """
    Test with BN-2d without treshold
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
    correction_indexes = get_fractionary_correction_indexes(vbm_df, treshold=0)
    assert correction_indexes["N"][0] == "p"
    assert correction_indexes["B"][0] == "p"


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
