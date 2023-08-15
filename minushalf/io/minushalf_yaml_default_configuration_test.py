"""
test minushalf_yaml_default_configuration module
"""
import numpy as np
from minushalf.io.minushalf_yaml_default_configuration import (
    VaspDefaultParams,
    CorrectionDefaultParams,
    AtomicProgramDefaultParams,
    MinushalfParams,
)


def test_to_list_atomic_program():
    """
    Test to_list method in atomic program params
    """
    params_list = AtomicProgramDefaultParams.to_list()
    assert params_list[0] == "pb"
    assert params_list[1] == "ae"
    assert params_list[2] == 100


def test_to_dict_atomic_program():
    """
    Test to_dict method in atomic program params
    """
    params_list = AtomicProgramDefaultParams.to_dict()
    assert params_list["exchange_correlation_code"] == "pb"
    assert params_list["calculation_code"] == "ae"
    assert params_list["max_iterations"] == 100


def test_to_list_vasp():
    """
    Test to_list method in vasp params
    """
    params_list = VaspDefaultParams.to_list()
    assert params_list[0] == ['mpirun', 'vasp']


def test_to_dict_vasp():
    """
    Test to_dict method in vasp params
    """
    params_list = VaspDefaultParams.to_dict()
    assert params_list["command"] == ['mpirun', 'vasp']


def test_to_list_correction():
    """
    test to_list method in correction params
    """
    params_list = CorrectionDefaultParams.to_list()
    assert len(params_list) == 11
    assert params_list[0] == "v"
    assert params_list[1] == "minushalf_potfiles"
    assert np.isclose(params_list[2], 1.0)
    assert params_list[3] is None
    assert params_list[4] is None
    assert np.isclose(params_list[5], 0.01)
    assert np.isclose(params_list[6], 10)
    assert np.isclose(params_list[7], 9)
    assert len(params_list[8]) == 0
    assert isinstance(params_list[8], list)
    assert len(params_list[9]) == 0
    assert isinstance(params_list[9], list)
    assert params_list[10] == False
    assert params_list[10] != None


def test_to_dict_correction():
    """
    test to_dict method in correction params
    """
    params_list = CorrectionDefaultParams.to_dict()
    assert params_list["inplace"] == False
    assert params_list["correction_code"] == "v"
    assert params_list["potfiles_folder"] == "minushalf_potfiles"
    assert np.isclose(params_list["amplitude"], 1.0)
    assert params_list["valence_cut_guess"] is None
    assert params_list["conduction_cut_guess"] is None
    assert np.isclose(params_list["tolerance"], 0.01)
    assert np.isclose(params_list["fractional_valence_treshold"], 10)
    assert np.isclose(params_list["fractional_conduction_treshold"], 9)
    assert isinstance(params_list["overwrite_vbm"], list)
    assert isinstance(params_list["overwrite_cbm"], list)


def test_to_list_minushalf_params():
    """
    test to_list method in minushalf params
    """
    params_list = MinushalfParams.to_list()
    assert params_list[0] == "software"
    assert params_list[1] == "atomic_program"
    assert params_list[2] == "correction"


def test_to_dict_minushalf_params():
    """
    test to_dict method in minushalf params
    """
    params_list = MinushalfParams.to_dict()
    assert params_list["software"] == "software"
    assert params_list["atomic_program"] == "atomic_program"
    assert params_list["correction"] == "correction"
