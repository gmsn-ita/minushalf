"""
Test software module
in data folder
"""
from minushalf.data import Softwares


def test_softwares():
    """
    Verify if all softwares available are
    in the class
    """
    softwares = ["VASP"]
    for element in softwares:
        assert Softwares[element.lower()].value == element


def test_softwares_get_default():
    """
    Method that returns the default value for this propertie
    """
    assert Softwares.get_default() == "VASP"
