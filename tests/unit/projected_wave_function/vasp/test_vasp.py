"""
Test file for vasp module for deal with
projeceted wave functions
"""
import pytest
import pandas as pd
from minushalf.projected_wave_function.vasp import Vasp


@pytest.mark.xfail(raises=FileNotFoundError)
def test_wrong_path_vasprun():
    """
    Test vasp creation wit
    """
    Vasp("./somewrongpath")


def test_vbm_charachter_with_gan(file_path):
    """
    Testing Vbm character function
    using the vasprun.xml relative to GaN.
    """
    expected_vbm_df = pd.read_csv(file_path("vbm_gan.csv"), index_col=0)

    vasp = Vasp(file_path("vasprun_GaN.xml"))
    vbm_character_df = vasp.character_vbm()

    assert vbm_character_df.equals(expected_vbm_df)


def test_cbm_character_with_gan(file_path):
    """
    Testing Cbm character function
    using the vasprun.xml relative to GaN.
    """
    expected_cbm_df = pd.read_csv(file_path("cbm_gan.csv"), index_col=0)

    vasp = Vasp(file_path("vasprun_GaN.xml"))
    cbm_character_df = vasp.character_cbm()

    assert cbm_character_df.equals(expected_cbm_df)
