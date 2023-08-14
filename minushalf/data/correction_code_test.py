"""
Test correction code module
in data folder
"""
from minushalf.data.correction_code import CorrectionCode


def test_correction_code():
    """
    Verify if all codes available are
    in the class
    """
    codes = ["v", "vf", "c", "cf", "vc", "vfc", "vcf", "vfcf"]
    for element in codes:
        assert CorrectionCode[element.lower()].value == element
