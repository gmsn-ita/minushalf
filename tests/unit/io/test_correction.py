"""
Test correction input class
"""
import numpy as np
from minushalf.io.correction import Correction


def test_default_parameters():
    """
    Test the default parameters 
    """
    correction = Correction()

    assert correction.correction_code == "v"
    assert correction.potfiles_folder == "minushalf_potfiles"
    assert np.isclose(correction.amplitude, 1.0)
    assert correction.valence_cut_guess == None
    assert correction.conduction_cut_guess == None
    assert np.isclose(correction.tolerance, 0.01)
    assert correction.fractional_valence_treshold == 10
    assert correction.fractional_conduction_treshold == 9
    assert correction.overwrite_cbm == []
    assert correction.overwrite_vbm == []
    assert correction.inplace == False


def test_override_parameters():
    """
    Test the class with other parameters
    """
    params = {
        "correction_code": "vf",
        "potfiles_folder": "potfiles",
        "amplitude": 2.0,
        "valence_cut_guess": [['C', 'p', 3.34]],
        "conduction_cut_guess": [['C', 'p', 3.34]],
        "tolerance": 0.1,
        "fractional_valence_treshold": 12,
        "fractional_conduction_treshold": 14,
        "overwrite_vbm": [1, 2],
        "overwrite_cbm": [1, 3],
        "inplace": True
    }
    correction = Correction(**params)
    assert correction.correction_code == "vf"
    assert correction.potfiles_folder == "potfiles"
    assert np.isclose(correction.amplitude, 2.0)
    assert correction.valence_cut_guess == [['C', 'p', 3.34]]
    assert correction.conduction_cut_guess == [['C', 'p', 3.34]]
    assert np.isclose(correction.tolerance, 0.1)
    assert correction.fractional_valence_treshold == 12
    assert correction.fractional_conduction_treshold == 14
    assert correction.overwrite_cbm == [1, 3]
    assert correction.overwrite_vbm == [1, 2]
    assert correction.inplace == True


def test_to_list():
    """
    Test method to_list
    """
    params_list = [
        "v", "minushalf_potfiles", 1.0, None, None, 0.01, 10, 9, [], [], None,
        None, False, None
    ]
    correction = Correction()

    assert correction.to_list() == params_list


def test_to_dict():
    """
    Test method to_list
    """
    params = {
        "correction_code": "vf",
        "potfiles_folder": "potfiles",
        "amplitude": 2.0,
        "valence_cut_guess": [['C', 'p', 3.34]],
        "conduction_cut_guess": [['C', 'p', 3.34]],
        "tolerance": 0.1,
        "fractional_valence_treshold": 12,
        "fractional_conduction_treshold": 14,
        "overwrite_vbm": [1, 2],
        "overwrite_cbm": [1, 3],
        "inplace": True,
        "vbm_characters": None,
        "cbm_characters": [["Ga", "d", "100"]],
        "divide_character": [["Ga", "p", 2]],
    }
    correction = Correction(**params)

    assert correction.to_dict() == params