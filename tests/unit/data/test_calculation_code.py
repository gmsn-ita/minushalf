"""
Test calculation code module
in data folder
"""
from minushalf.data.calculation_code import CalculationCode


def test_calculation_code():
    """
    Verify if all codes available are
    in the class
    """
    calculation_code = ["ae"]
    for element in calculation_code:
        assert CalculationCode[element.lower()].value == element


def test_calculation_code_get_default():
    """
    Method that returns the default value for this propertie
    """
    assert CalculationCode.get_default() == "ae"
