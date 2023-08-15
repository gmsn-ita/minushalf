"""
Test run-atomic command in minushalf CLI
"""
import os
from click.testing import CliRunner
from minushalf import atomic_program
from minushalf.commands.create_input import create_input
from minushalf.utils.electronic_distribution import ElectronicDistribution
from minushalf.io.input_file import InputFile


def test_create_input_with_all_elements():
    """
    Test if the atomic run produce the correct output for
    elements in the periodic table
    """
    runner = CliRunner()

    for element in ElectronicDistribution:
        symbol = str(element).upper()
        with runner.isolated_filesystem():
            result = runner.invoke(create_input, [symbol])
            assert result.exit_code == 0
            assert os.path.exists("INP") == True
