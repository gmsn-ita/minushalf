"""
Test occupation command in minushalf CLI
"""
import numpy as np
import pytest
from click.testing import CliRunner
from minushalf import atomic_program
from minushalf.io import InputFile
from minushalf.commands import occupation


def test_occupation_with_default_params_O():
    """
    Test occupation command in the p orbital of oxygen
    """
    inp_oxygen = InputFile.minimum_setup("O", "pb")
    lines = inp_oxygen.to_stringlist()
    runner = CliRunner()
    with runner.isolated_filesystem():
        with open("INP", "w") as file:
            file.writelines(lines)

        result = runner.invoke(occupation, ["1"])
        assert result.exit_code == 0
        occupied_file = InputFile.from_file("INP_OCC")
        for orbital in occupied_file.valence_orbitals:
            if orbital["l"] == 1:
                assert np.isclose(orbital["occupation"][0], 3.5) == True


def test_multiple_occupation_with_default_params_O():
    """
    Test occupation command in the p and s orbitals of oxygen
    """
    inp_oxygen = InputFile.minimum_setup("O", "pb")
    lines = inp_oxygen.to_stringlist()
    runner = CliRunner()
    with runner.isolated_filesystem():
        with open("INP", "w") as file:
            file.writelines(lines)

        result = runner.invoke(occupation, ["0,1", "*,*"])
        assert result.exit_code == 0
        occupied_file = InputFile.from_file("INP_OCC")
        for orbital in occupied_file.valence_orbitals:
            if orbital["l"] == 1:
                assert np.isclose(orbital["occupation"][0], 3.5) == True
            if orbital["l"] == 0:
                assert np.isclose(orbital["occupation"][0], 1.5) == True


def test_failing_occupation_with_default_params_O():
    """
    Test passing invalid secondary quantum number
    """
    inp_oxygen = InputFile.minimum_setup("O", "pb")
    lines = inp_oxygen.to_stringlist()
    runner = CliRunner()
    with runner.isolated_filesystem():
        with open("INP", "w") as file:
            file.writelines(lines)

        result = runner.invoke(occupation, ["0,-1", "*,*"])
        assert result.exit_code == 1


def test_occupation_with_default_params_Yb():
    """
    Test occupation command in the p orbital of ytterbium
    """
    inp_oxygen = InputFile.minimum_setup("Yb", "pb")
    lines = inp_oxygen.to_stringlist()
    runner = CliRunner()
    with runner.isolated_filesystem():
        with open("INP", "w") as file:
            file.writelines(lines)

        result = runner.invoke(occupation, ["0", "50"])
        assert result.exit_code == 0
        occupied_file = InputFile.from_file("INP_OCC")
        for orbital in occupied_file.valence_orbitals:
            if orbital["l"] == 0:
                assert np.isclose(orbital["occupation"][0], 1.75) == True


def test_multiple_occupation_with_default_params_Yb():
    """
    Test multiple occupation command in the p orbital of ytterbium
    """
    inp_oxygen = InputFile.minimum_setup("Yb", "pb")
    lines = inp_oxygen.to_stringlist()
    runner = CliRunner()
    with runner.isolated_filesystem():
        with open("INP", "w") as file:
            file.writelines(lines)

        result = runner.invoke(occupation, ["0,3", "50, 50"])
        assert result.exit_code == 0
        occupied_file = InputFile.from_file("INP_OCC")
        for orbital in occupied_file.valence_orbitals:
            if orbital["l"] == 0:
                assert np.isclose(orbital["occupation"][0], 1.75) == True
            if orbital["l"] == 3:
                assert np.isclose(orbital["occupation"][0], 13.75) == True


def test_failing_occupation_with_default_params_Yb():
    """
    Test invalid occupation number for Yb
    """
    inp_oxygen = InputFile.minimum_setup("Yb", "pb")
    lines = inp_oxygen.to_stringlist()
    runner = CliRunner()
    with runner.isolated_filesystem():
        with open("INP", "w") as file:
            file.writelines(lines)

        result = runner.invoke(occupation, ["0,3", "50,101"])
        assert result.exit_code == 1


def test_occupation_with_default_params_Na():
    """
    Test occupation command in the p orbital of sodium
    """
    inp_oxygen = InputFile.minimum_setup("Na", "pb")
    lines = inp_oxygen.to_stringlist()
    runner = CliRunner()
    with runner.isolated_filesystem():
        with open("INP", "w") as file:
            file.writelines(lines)

        result = runner.invoke(occupation, ["0", "80"])
        assert result.exit_code == 0
        occupied_file = InputFile.from_file("INP_OCC")
        for orbital in occupied_file.valence_orbitals:
            if orbital["l"] == 0:
                assert np.isclose(orbital["occupation"][0], 0.6) == True


def test_failing_occupation_with_default_params_Na():
    """
    Test wrong input to secondary quantum number
    """
    inp_oxygen = InputFile.minimum_setup("Na", "pb")
    lines = inp_oxygen.to_stringlist()
    runner = CliRunner()
    with runner.isolated_filesystem():
        with open("INP", "w") as file:
            file.writelines(lines)

        result = runner.invoke(occupation, ["4", "80"])
        assert result.exit_code == 1


def test_occupation_with_default_params_Au():
    """
    Test occupation command in the p orbital of gold
    """
    inp_oxygen = InputFile.minimum_setup("Au", "pb")
    lines = inp_oxygen.to_stringlist()
    runner = CliRunner()
    with runner.isolated_filesystem():
        with open("INP", "w") as file:
            file.writelines(lines)

        result = runner.invoke(occupation, ["2", "100"])
        assert result.exit_code == 0
        occupied_file = InputFile.from_file("INP_OCC")
        for orbital in occupied_file.valence_orbitals:
            if orbital["l"] == 2:
                assert np.isclose(orbital["occupation"][0], 9.5) == True


def test_failing_occupation_with_default_params_Au():
    """
    Test wrong occuparion percentual passed to command
    """
    inp_oxygen = InputFile.minimum_setup("Au", "pb")
    lines = inp_oxygen.to_stringlist()
    runner = CliRunner()
    with runner.isolated_filesystem():
        with open("INP", "w") as file:
            file.writelines(lines)

        result = runner.invoke(occupation, ["2", "-1"])
        assert result.exit_code == 2


def test_occupation_with_default_params_Si():
    """
    Test occupation command in the p orbital of silicium
    """
    inp_oxygen = InputFile.minimum_setup("Si", "pb")
    lines = inp_oxygen.to_stringlist()
    runner = CliRunner()
    with runner.isolated_filesystem():
        with open("INP", "w") as file:
            file.writelines(lines)

        result = runner.invoke(occupation, ["1", "0"])
        assert result.exit_code == 0
        occupied_file = InputFile.from_file("INP_OCC")
        for orbital in occupied_file.valence_orbitals:
            if orbital["l"] == 1:
                assert np.isclose(orbital["occupation"][0], 2.0) == True
