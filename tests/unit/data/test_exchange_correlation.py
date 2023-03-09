"""
Test exchange correlation module
in data folder
"""
from minushalf.data.exchange_correlation import ExchangeCorrelation


def test_exchange_correlation():
    """
    Verify if all codes available are
    in the class
    """
    exchange_code = ["ca", "wi", "hl", "gl", "bh", "pb", "rp", "rv", "bl"]
    for element in exchange_code:
        assert ExchangeCorrelation[element.lower()].value == element


def test_exchange_correlation_get_default():
    """
    Method that returns the default value for this propertie
    """
    assert ExchangeCorrelation.get_default() == "pb"