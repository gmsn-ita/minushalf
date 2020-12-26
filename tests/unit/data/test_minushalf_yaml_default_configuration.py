"""
test minushalf_yaml_default_configuration module
"""
from minushalf.data import (VaspDefaultParams, CorrectionDefaultParams,
                            AtomicProgramDefaultParams)


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


def test_to_dict_correction():
    """
    test to_dict method in correction params
    """
    params_list = CorrectionDefaultParams.to_dict()
    assert params_list["correction_code"] == "v"
