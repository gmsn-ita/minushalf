"""
Test proj_output module
The functions in this file obey the following rules:
test_(what the function is meant to do)_(chemical compound)
"""
import numpy as np
from minushalf.softwares.qe.projoutput import ProjOutput


def test_parse_projoutput_header_aln_2d(file_path):
    """
    Test if the dimensions and states_info are correctly parsed
    from the projwfc_up file for AlN (2D).
    """
    filename = file_path("aln-2d/proj.projwfc_up")
    proj = ProjOutput(filename)

    # --- dimensions ---
    assert proj.num_kpoints == 40   
    assert proj.num_bands   == 20   
    assert proj.num_states  == 16   

    # --- states_info: first state ---
    first = proj.states_info[0]
    assert first["state_idx"] == 1
    assert first["atom_idx"]  == 1   
    assert first["symbol"]    == "Al"
    assert first["l"]         == 0
    assert first["m"]         == 1

    # --- states_info: last state ---
    last = proj.states_info[-1]
    assert last["state_idx"] == proj.num_states
    assert last["atom_idx"] == 4   
    assert last["symbol"]    == "N"
    assert last["l"]        == 1
    assert last["m"]        == 3



def test_get_band_projection_kpt_1_band_1_aln_2d(file_path):
    """
    Verify the band projection for kpoint 1, band 1 for AlN (2D).
    """
    filename = file_path("aln-2d/proj.projwfc_up")
    proj = ProjOutput(filename)

    # Expected: { atom_index_str: [orbital_proj_0, orbital_proj_1, ...] }
    expected = {
         "1": [0.4379301244, 0.000186398, 0.0, 0.0],
         "2": [0.4379301244, 0.000186398, 0.0, 0.0],
         "3": [0.0012675996, 3.4402e-06, 0.0, 0.0],
         "4": [0.0012675996, 3.4402e-06, 0.0, 0.0],
    }

    band_projection = proj.get_band_projection(kpoint=1, band_number=1)

    # Check that the returned keys match exactly
    assert set(band_projection.keys()) == set(expected.keys())

    for atom_index, projections in band_projection.items():
        # Check that this atom_index was expected
        assert atom_index in expected

        # Check that the number of orbitals matches
        assert len(projections) == len(expected[atom_index])

        # Check each orbital value
        for orbital_idx, value in enumerate(projections):
            assert np.isclose(value, expected[atom_index][orbital_idx]), (
                f"atom {atom_index}, orbital {orbital_idx}: "
                f"got {value}, expected {expected[atom_index][orbital_idx]}"
            )


def test_get_band_projection_kpt_1_band_last_aln_2d(file_path):
    """
    Verify the band projection for kpoint 1, last band for AlN (2D).
    """
    filename = file_path("aln-2d/proj.projwfc_up")
    proj = ProjOutput(filename)

    expected = {
         "1": [0.0, 0.0, 0.0005205484, 0.0005205484],
         "2": [0.0, 0.0, 0.0005205484, 0.0005205484],
         "3": [0.0, 0.0, 0.0011697102, 0.0011697102],
         "4": [0.0, 0.0, 0.0011697102, 0.0011697102],
    }

    band_projection = proj.get_band_projection(
        kpoint=1, band_number=proj.num_bands
    )

    # Check that the returned keys match exactly
    assert set(band_projection.keys()) == set(expected.keys())

    for atom_index, projections in band_projection.items():
        # Check that this atom_index was expected
        assert atom_index in expected

        # Check that the number of orbitals matches
        assert len(projections) == len(expected[atom_index])

        # Check each orbital value
        for orbital_idx, value in enumerate(projections):
            assert np.isclose(value, expected[atom_index][orbital_idx]), (
                f"atom {atom_index}, orbital {orbital_idx}: "
                f"got {value}, expected {expected[atom_index][orbital_idx]}"
            )


def test_get_band_projection_kpt_last_band_1_aln_2d(file_path):
    """
    Verify the band projection for the last kpoint, band 1 for AlN (2D).
    """
    filename = file_path("aln-2d/proj.projwfc_up")
    proj = ProjOutput(filename)

    expected = {
         "1": [0.0335613473, 0.0695273437, 0.1269512691, 0.1269512691],
         "2": [0.0335613473, 0.0695273437, 0.1269512691, 0.1269512691],
         "3": [0.0006193388, 0.0009559325, 0.0007989039, 0.0007989039],
         "4": [0.0006193388, 0.0009559325, 0.0007989039, 0.0007989039],
    }

    band_projection = proj.get_band_projection(
        kpoint=proj.num_kpoints, band_number=1
    )

    # Check that the returned keys match exactly
    assert set(band_projection.keys()) == set(expected.keys())

    for atom_index, projections in band_projection.items():
        # Check that this atom_index was expected
        assert atom_index in expected

        # Check that the number of orbitals matches
        assert len(projections) == len(expected[atom_index])

        # Check each orbital value
        for orbital_idx, value in enumerate(projections):
            assert np.isclose(value, expected[atom_index][orbital_idx]), (
                f"atom {atom_index}, orbital {orbital_idx}: "
                f"got {value}, expected {expected[atom_index][orbital_idx]}"
            )

def test_parse_projoutput_header_sic_2d(file_path):
    """
    Test if the dimensions and states_info are correctly parsed
    from the projwfc_up file for SiC (2D).
    """
    filename = file_path("sic-2d/proj.projwfc_up")
    proj = ProjOutput(filename)

    # --- dimensions ---
    assert proj.num_kpoints == 29   
    assert proj.num_bands   == 8   
    assert proj.num_states  == 8   

    # --- states_info: first state ---
    first = proj.states_info[0]
    assert first["state_idx"] == 1
    assert first["atom_idx"]  == 1   
    assert first["symbol"]    == "Si"
    assert first["l"]         == 0
    assert first["m"]         == 1

    # --- states_info: last state ---
    last = proj.states_info[-1]
    assert last["state_idx"] == proj.num_states
    assert last["atom_idx"] == 2   
    assert last["symbol"]    == "C"
    assert last["l"]        == 1
    assert last["m"]        == 3



def test_get_band_projection_kpt_1_band_1_sic_2d(file_path):
    """
    Verify the band projection for kpoint 1, band 1 for SiC (2D).
    """
    filename = file_path("sic-2d/proj.projwfc_up")
    proj = ProjOutput(filename)

    # Expected: { atom_index_str: [orbital_proj_0, orbital_proj_1, ...] }
    expected = {
         "1": [0.4879185383, 0.0, 0.0, 0.0],
         "2": [0.5099774512, 0.0, 0.0, 0.0],
    }

    band_projection = proj.get_band_projection(kpoint=1, band_number=1)

    # Check that the returned keys match exactly
    assert set(band_projection.keys()) == set(expected.keys())

    for atom_index, projections in band_projection.items():
        # Check that this atom_index was expected
        assert atom_index in expected

        # Check that the number of orbitals matches
        assert len(projections) == len(expected[atom_index])

        # Check each orbital value
        for orbital_idx, value in enumerate(projections):
            assert np.isclose(value, expected[atom_index][orbital_idx]), (
                f"atom {atom_index}, orbital {orbital_idx}: "
                f"got {value}, expected {expected[atom_index][orbital_idx]}"
            )


def test_get_band_projection_kpt_1_band_last_sic_2d(file_path):
    """
    Verify the band projection for kpoint 1, last band for SiC (2D).
    """
    filename = file_path("sic-2d/proj.projwfc_up")
    proj = ProjOutput(filename)

    expected = {
         "1": [0.0, 0.2882750648, 0.2882750648, 0.2882750648],
         "2": [0.0, 0.0376281208, 0.0376281208, 0.0376281208],
    }

    band_projection = proj.get_band_projection(
        kpoint=1, band_number=proj.num_bands
    )

    # Check that the returned keys match exactly
    assert set(band_projection.keys()) == set(expected.keys())

    for atom_index, projections in band_projection.items():
        # Check that this atom_index was expected
        assert atom_index in expected

        # Check that the number of orbitals matches
        assert len(projections) == len(expected[atom_index])

        # Check each orbital value
        for orbital_idx, value in enumerate(projections):
            assert np.isclose(value, expected[atom_index][orbital_idx]), (
                f"atom {atom_index}, orbital {orbital_idx}: "
                f"got {value}, expected {expected[atom_index][orbital_idx]}"
            )


def test_get_band_projection_kpt_last_band_1_sic_2d(file_path):
    """
    Verify the band projection for the last kpoint, band 1 for SiC (2D).
    """
    filename = file_path("sic-2d/proj.projwfc_up")
    proj = ProjOutput(filename)

    expected = {
         "1": [0.0, 0.0906329462, 0.0906329462, 0.0906329462],
         "2": [0.7243652644, 0.0, 0.0, 0.0],
    }

    band_projection = proj.get_band_projection(
        kpoint=proj.num_kpoints, band_number=1
    )

    # Check that the returned keys match exactly
    assert set(band_projection.keys()) == set(expected.keys())

    for atom_index, projections in band_projection.items():
        # Check that this atom_index was expected
        assert atom_index in expected

        # Check that the number of orbitals matches
        assert len(projections) == len(expected[atom_index])

        # Check each orbital value
        for orbital_idx, value in enumerate(projections):
            assert np.isclose(value, expected[atom_index][orbital_idx]), (
                f"atom {atom_index}, orbital {orbital_idx}: "
                f"got {value}, expected {expected[atom_index][orbital_idx]}"
            )
