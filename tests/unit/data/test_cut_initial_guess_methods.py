"""
Test the methods available to give cut initial gesses
"""
from minushalf.data import CutInitialGuessMethods


def test_cut_initial_guess():
    """
    Verify if all methods available are
    in the class
    """
    methods = ["three_dimensions"]
    values = ["3d"]
    for element, value in zip(methods, values):
        assert CutInitialGuessMethods[element.lower()].value == value


def test_to_list_method():
    """
    Test method to list
    """
    methods = CutInitialGuessMethods.to_list()
    assert len(methods) == 1
    assert methods[0] == "3d"
