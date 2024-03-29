"""
Test run-atomic command in minushalf CLI
"""
import os
from unittest import mock
from click.testing import CliRunner
from minushalf import atomic_program
from minushalf.commands.run_atomic_program import run_atomic
from minushalf.utils.electronic_distribution import ElectronicDistribution
from minushalf.io.input_file import InputFile


def test_atomic_run_workflow():
    """
    Test the workflow of the command run-atomic
    """
    with mock.patch("minushalf.atomic_program.run"):
        runner = CliRunner()
        runner.invoke(run_atomic, [])
        atomic_program.run.assert_called_once()


def test_atomic_with_all_elements():
    """
    Test if the atomic run produce the correct output for
    elements in the periodic table
    """
    runner = CliRunner()

    for element in ElectronicDistribution:
        inp = InputFile.minimum_setup(str(element), "pb")
        lines = inp.to_stringlist()
        with runner.isolated_filesystem():
            with open("INP", "w") as file:
                file.writelines(lines)
            result = runner.invoke(run_atomic, [])
            assert result.exit_code == 0
            assert os.path.exists("VTOTAL0") == True
            assert os.path.exists("VTOTAL.ae") == True
            assert os.path.exists("VTOTAL2") == True
            assert os.path.exists("VTOTAL3") == True
