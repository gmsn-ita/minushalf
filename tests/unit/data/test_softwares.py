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
