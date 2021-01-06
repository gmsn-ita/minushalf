"""
test minushalf_yaml_default_configuration module
"""
import numpy as np
from minushalf.data import (
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
    assert params_list[0] == 4
    assert params_list[1] == "vasp"


def test_to_dict_vasp():
    """
    Test to_dict method in vasp params
    """
    params_list = VaspDefaultParams.to_dict()
    assert params_list["number_of_cores"] == 4
    assert params_list["path"] == "vasp"


def test_to_list_correction():
    """
    test to_list method in correction params
    """
    params_list = CorrectionDefaultParams.to_list()
    assert params_list[0] == "v"
    assert params_list[1] == "minushalf_potfiles"
    assert np.isclose(params_list[2], 1.0)


def test_to_dict_correction():
    """
    test to_dict method in correction params
    """
    params_list = CorrectionDefaultParams.to_dict()
    assert params_list["correction_code"] == "v"
    assert params_list["potfiles_folder"] == "minushalf_potfiles"
    assert np.isclose(params_list["amplitude"], 1.0)


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
