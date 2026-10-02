"""
Test create-input command in minushalf CLI
"""
import os
from click.testing import CliRunner
from minushalf.commands.create_input import create_input
from minushalf.utils.electronic_distribution import ElectronicDistribution


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

def test_create_input_with_all_elements_qe():
    """
    Test if the atomic run produce the correct output for
    elements in the periodic table
    """
    runner = CliRunner()

    for element in ElectronicDistribution:
        symbol = str(element).upper()
        with runner.isolated_filesystem():
            result = runner.invoke(create_input, [symbol, '-s', 'QE'])
            assert result.exit_code == 0
            assert os.path.exists("INP") == True

