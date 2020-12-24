"""
Test run-atomic command in minushalf CLI
"""
import os
from click.testing import CliRunner
from minushalf import atomic_program
from minushalf.commands import create_input
from minushalf.data import ElectronicDistribution
from minushalf.utils import InputFile


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
