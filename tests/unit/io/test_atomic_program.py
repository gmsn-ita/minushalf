"""
Test atomic program input class
"""
from minushalf.io import AtomicProgram


def test_default_parameters():
    """
    Test the default parameters 
    """
    atomic_program = AtomicProgram()

    assert atomic_program.exchange_correlation_code == "pb"
    assert atomic_program.calculation_code == "ae"
    assert atomic_program.max_iterations == 100


def test_override_parameters():
    """
    Test the class with other parameters
    """
    atomic_program = AtomicProgram(exchange_correlation_code="wi",
                                   calculation_code="he",
                                   max_iterations=300)

    assert atomic_program.exchange_correlation_code == "wi"
    assert atomic_program.calculation_code == "he"
    assert atomic_program.max_iterations == 300


def test_to_list():
    """
    Test method to_list
    """
    atomic_program = AtomicProgram()

    parameters_list = atomic_program.to_list()

    assert parameters_list == ["pb", "ae", 100]


def test_to_dict():
    """
    Test method to_list
    """
    atomic_program = AtomicProgram()

    parameters_dict = atomic_program.to_dict()

    assert parameters_dict == {
        "exchange_correlation_code": "pb",
        "calculation_code": "ae",
        "max_iterations": 100
    }
