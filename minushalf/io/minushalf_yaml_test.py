"""
Test the class that leads with the minushalf.yaml file,
the input file for the execute command
"""
import pytest
from minushalf.io.minushalf_yaml import MinushalfYaml
from minushalf.softwares.softwares import (Softwares)


def test_default_parameters():
    """
    Test minushalf.yaml class withou pass
    parameters, so all parameters
    are set to default
    """
    file = MinushalfYaml.from_file()
    params = {
        **file.get_atomic_program_params(),
        **file.get_correction_params(),
        **file.get_software_configurations_params()
    }
    default_params = {
        "exchange_correlation_code": "pb",
        "calculation_code": "ae",
        "max_iterations": 100,
        "potfiles_folder": "minushalf_potfiles",
        "amplitude": 1.0,
        "valence_cut_guess": None,
        "conduction_cut_guess": None,
        "tolerance": 0.01,
        "fractional_valence_treshold": 10,
        "fractional_conduction_treshold": 9,
        "indirect": False,
        "correction_code": "v",
        "overwrite_vbm": [],
        "overwrite_cbm": [],
        "command": ["mpirun", "vasp"],
        "vbm_characters": None,
        "cbm_characters": None,
        "divide_character": None
    }
    assert file.get_command() == ["mpirun", "vasp"]
    assert file.get_correction_code() == "v"
    assert file.get_overwrite_cbm() == []
    assert file.get_overwrite_vbm() == []
    assert params == default_params


def test_minushalf_without_filling_correction(file_path):
    """
    Test minushalf.yaml class without filling
    correction section in minushalf.yaml
    """
    minushalf_path = file_path(
        "/minushalf_yaml/minushalf_partially_filled.yaml")
    file = MinushalfYaml.from_file(minushalf_path)
    params = {
        **file.get_atomic_program_params(),
        **file.get_correction_params(),
        **file.get_software_configurations_params()
    }
    expected_params = {
        "exchange_correlation_code": "wi",
        "calculation_code": "ae",
        "max_iterations": 200,
        "potfiles_folder": "minushalf_potfiles",
        "amplitude": 1.0,
        "valence_cut_guess": None,
        "conduction_cut_guess": None,
        "tolerance": 0.01,
        "fractional_valence_treshold": 10,
        "fractional_conduction_treshold": 9,
        "indirect": False,
        "correction_code": "v",
        "overwrite_vbm": [],
        "overwrite_cbm": [],
        "command": ["mpirun", "-np", "6", "../vasp"],
        "vbm_characters": None,
        "cbm_characters": None,
        "divide_character": None
    }
    assert file.get_command() == ["mpirun", "-np", "6", "../vasp"]
    assert file.get_software_name() == Softwares.vasp.value
    assert file.get_correction_code() == "v"
    assert file.get_overwrite_cbm() == []
    assert file.get_overwrite_vbm() == []
    assert params == expected_params


def test_minushalf_filled_out(file_path):
    """
    Test minushalf.yaml class with all
    parameters modified
    """
    minushalf_path = file_path("/minushalf_yaml/minushalf_filled_out.yaml")
    file = MinushalfYaml.from_file(minushalf_path)
    params = {
        **file.get_atomic_program_params(),
        **file.get_correction_params(),
        **file.get_software_configurations_params()
    }
    expected_params = {
        "exchange_correlation_code": "wi",
        "calculation_code": "ae",
        "max_iterations": 200,
        "potfiles_folder": "../potcar",
        "amplitude": 3.0,
        "valence_cut_guess": [["C", "p", 3.45]],
        "conduction_cut_guess": [["C", "p", 1.0]],
        "tolerance": 0.001,
        "fractional_valence_treshold": 15,
        "fractional_conduction_treshold": 23,
        "indirect": True,
        "correction_code": "vf",
        "overwrite_vbm": [1, 3],
        "overwrite_cbm": [1, 4],
        "vbm_characters": None,
        "cbm_characters": None,
        "command": ["mpirun", "-np", "6", "../vasp"],
        "divide_character": None
    }
    assert file.get_command() == ["mpirun", "-np", "6", "../vasp"]
    assert file.get_software_name() == Softwares.vasp.value
    assert file.get_correction_code() == "vf"
    assert file.get_overwrite_cbm() == [1, 4]
    assert file.get_overwrite_vbm() == [1, 3]
    assert file.get_vbm_characters() == None
    assert file.get_cbm_characters() == None
    assert params == expected_params


def test_minushalf_characters(file_path):
    """
    Test minushalf.yaml class with all
    parameters modified
    """
    minushalf_path = file_path("/minushalf_yaml/minushalf_characters.yaml")
    file = MinushalfYaml.from_file(minushalf_path)
    params = {
        **file.get_atomic_program_params(),
        **file.get_correction_params(),
        **file.get_software_configurations_params()
    }
    expected_params = {
        "exchange_correlation_code": "wi",
        "calculation_code": "ae",
        "max_iterations": 200,
        "potfiles_folder": "../potcar",
        "amplitude": 3.0,
        "valence_cut_guess": [["C", "p", 3.45]],
        "conduction_cut_guess": [["C", "p", 1.0]],
        "tolerance": 0.001,
        "fractional_valence_treshold": 15,
        "fractional_conduction_treshold": 23,
        "indirect": True,
        "correction_code": "vf",
        "overwrite_vbm": [1, 3],
        "overwrite_cbm": [1, 4],
        "vbm_characters": [["C", "p", 23.0]],
        "cbm_characters": [["Ga", "d", 23.0]],
        "command": ["mpirun", "-np", "6", "../vasp"],
        "divide_character": None
    }
    assert file.get_command() == ["mpirun", "-np", "6", "../vasp"]
    assert file.get_software_name() == Softwares.vasp.value
    assert file.get_correction_code() == "vf"
    assert file.get_overwrite_cbm() == [1, 4]
    assert file.get_overwrite_vbm() == [1, 3]
    assert file.get_vbm_characters() == [["C", "p", 23.0]]
    assert file.get_cbm_characters() == [["Ga", "d", 23.0]]
    assert file.get_potential_folder() == "../potcar"
    assert file.get_amplitude() == 3.0
    assert file.get_max_iterations() == 200
    assert file.get_exchange_corr_code() == "wi"
    assert file.get_calculation_code() == "ae"
    assert file.get_valence_cut_initial_guess() == [["C", "p", 3.45]]
    assert file.get_conduction_cut_initial_guess() == [["C", "p", 1.0]]
    assert file.get_tolerance() == 0.001
    assert file.get_indirect() == True
    assert file.get_divide_character() == None
    assert params == expected_params


@pytest.mark.xfail
def test_minushalf_wrong_software_name(file_path):
    """
    Test minushalf.yaml passing
    a wrong software name
    """
    minushalf_path = file_path(
        "/minushalf_yaml/minushalf_invalid_software.yaml")
    MinushalfYaml.from_file(minushalf_path)


@pytest.mark.xfail
def test_minushalf_wrong_correction_code_name(file_path):
    """
    Test minushalf.yaml passing
    a invalid correction code
    """
    minushalf_path = file_path(
        "/minushalf_yaml/minushalf_invalid_correction_code.yaml")
    MinushalfYaml.from_file(minushalf_path)
