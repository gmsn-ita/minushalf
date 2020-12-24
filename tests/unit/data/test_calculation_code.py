"""
Test calculation code module
in data folder
"""
from minushalf.data import CalculationCode


def test_calculation_code():
    """
    Verify if all codes available are
    in the class
    """
    calculation_code = ["ae"]
    for element in calculation_code:
        assert CalculationCode[element.lower()].value == element
