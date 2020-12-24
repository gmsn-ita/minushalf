"""
Test exchange correlation module
in data folder
"""
from minushalf.data import ExchangeCorreltion


def test_exchange_correlation():
    """
    Verify if all codes available are
    in the class
    """
    exchange_code = ["ca", "wi", "hl", "gl", "bh", "pb", "rp", "rv", "bl"]
    for element in exchange_code:
        assert ExchangeCorreltion[element.lower()].value == element
