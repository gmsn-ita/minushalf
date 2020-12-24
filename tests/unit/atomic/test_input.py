"""
Test input file class
"""
import pytest
import numpy as np
from collections import Counter
from minushalf.utils import ElectronicDistribution, PeriodicTable
from minushalf.atomic import InputFile


def test_minimum_setup(file_path):
    """
    Test minimum setup function to
    generate INP files with the elements
    symbol and the functional of exchange and exchange
    correlation
    """

    for element in ElectronicDistribution:
        symbol = str(element)
        path = file_path(f"{symbol}/INP")
        inp = InputFile.minimum_setup(symbol, 'pb')
        with open(path, "r") as file:
            assert file.read() == "".join(inp.to_stringlist())


def test_from_file(file_path):
    """
    Test from_file path to generate instances
    of the InputFile class with an INP text file
    """
    for element in ElectronicDistribution:
        symbol = str(element)
        path = file_path(f"{symbol}/INP_COMMENTED")

        inp_from_file = InputFile.from_file(path)
        inp_minimum_setup = InputFile.minimum_setup(symbol, 'pb')

        from_file_string = "".join(inp_from_file.to_stringlist())
        minimum_setup_string = "".join(inp_minimum_setup.to_stringlist())

        assert from_file_string == minimum_setup_string


def test_electron_occupation_Ag():
    """
        Test occupation in INP file for Ag
        """
    inp = InputFile.minimum_setup('Ag', 'pb')
    inp.electron_occupation(0.5, 2)
    inp.electron_occupation(0.2, 0)

    for orbital in inp.valence_orbitals:
        if orbital["l"] == 2:
            assert np.isclose(orbital["occupation"][0], 9.5) == True
        elif orbital["l"] == 0:
            assert np.isclose(orbital["occupation"][0], 0.8) == True


def test_electron_occupation_Si():
    """
        Test occupation in INP file for Si
        """
    inp = InputFile.minimum_setup('Si', 'pb')
    inp.electron_occupation(0.4, 0)
    inp.electron_occupation(0.3, 1)

    for orbital in inp.valence_orbitals:
        if orbital["l"] == 0:
            assert np.isclose(orbital["occupation"][0], 1.6) == True
        elif orbital["l"] == 1:
            assert np.isclose(orbital["occupation"][0], 1.7) == True


def test_electron_occupation_Au():
    """
        Test occupation in INP file for Au
        """
    inp = InputFile.minimum_setup('Au', 'pb')
    inp.electron_occupation(0.1, 0)
    inp.electron_occupation(0.4, 2)

    for orbital in inp.valence_orbitals:
        if orbital["l"] == 0:
            assert np.isclose(orbital["occupation"][0], 0.9) == True
        elif orbital["l"] == 2:
            assert np.isclose(orbital["occupation"][0], 9.6) == True


def test_electron_occupation_Yb():
    """
        Test occupation in INP file for Yb
        """
    inp = InputFile.minimum_setup('Yb', 'pb')
    inp.electron_occupation(0.5, 0)
    inp.electron_occupation(0.5, 3)

    for orbital in inp.valence_orbitals:
        if orbital["l"] == 0:
            assert np.isclose(orbital["occupation"][0], 1.5) == True
        elif orbital["l"] == 3:
            assert np.isclose(orbital["occupation"][0], 13.5) == True


def test_electron_occupation_Na():
    """
        Test occupation in INP file for Na
        """
    inp = InputFile.minimum_setup('Na', 'pb')
    inp.electron_occupation(0.5, 0)

    for orbital in inp.valence_orbitals:
        if orbital["l"] == 0:
            assert np.isclose(orbital["occupation"][0], 0.5) == True


def test_exchange_and_correlation_functional():
    """
        Test if the factor of correlation and exchange
        will raises exception in the class construction.
        """
    exchange_correlation_types = [
        "ca", "wi", "hl", "gl", "bh", "pb", "rp", "rv", "bl"
    ]
    for type_correlation in exchange_correlation_types:
        InputFile(type_correlation, "ae", "Ge", "", 1, 1, [])


@pytest.mark.xfail
def test_exchange_and_correlation_functional_relativistic():
    """
        Test if the factor of correlation and exchange with
        relativistic factor will raises exception in the class construction.
        """
    exchange_correlation_types = [
        "car", "wir", "hlr", "glr", "bhr", "pbr", "rpr", "rvr", "blr"
    ]
    for type_correlation in exchange_correlation_types:
        InputFile(type_correlation, "ae", "Ge", "", 1, 1, [])


@pytest.mark.xfail
def test_exchange_and_correlation_functional_spin():
    """
        Test if the factor of correlation and exchange with
        spin factor will raises exception in the class construction.
        """
    exchange_correlation_types = [
        "cas", "wis", "hls", "gls", "bhs", "pbs", "rps", "rvs", "bls"
    ]
    for type_correlation in exchange_correlation_types:
        InputFile(type_correlation, "ae", "Ge", "", 1, 1, [])


def test_chemical_symbol():
    """
        Test if the chemical symbol of the predetermined
        will raises exception in the class construction.
        """
    for element in PeriodicTable:
        InputFile("pb", "ae", element.value, "", 1, 1, [])


@pytest.mark.xfail
def test_pass_wrong_exchange_functional():
    """
        Pass wrong correlation and exchange functional and
        expect to fail
        """
    InputFile('sd', 'ae', 'Ge', "", 1, 1, [])


@pytest.mark.xfail
def test_pass_wrong_calculation_code():
    """
        Pass wrong calculation code and
        expect to fail
        """
    InputFile('pb', 'ee', 'Ge', "", 1, 1, [])


@pytest.mark.xfail
def test_pass_wrong_chemical_symbol():
    """
        Pass wrong chemical symbol and
        expect to fail
        """
    InputFile('pb', 'ae', 'Ss', "", 1, 1, [])
