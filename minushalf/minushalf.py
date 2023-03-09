"""
Definition of the minushalf CLI
"""
import click
from minushalf.commands.run_atomic_program import run_atomic
from minushalf.commands.fractional_occupation import occupation
from minushalf.commands.vbm_character import vbm_character
from minushalf.commands.cbm_character import cbm_character
from minushalf.commands.band_character import band_character
from minushalf.commands.create_input import create_input
from minushalf.commands.correct_potfile import correct_potfile
from minushalf.commands.band_gap import band_gap
from minushalf.commands.execute import execute


@click.group()
def minushalf():
    """
    CLI for automating DFT -1/2 calculations
    """


minushalf.add_command(run_atomic)
minushalf.add_command(occupation)
minushalf.add_command(vbm_character)
minushalf.add_command(cbm_character)
minushalf.add_command(band_character)
minushalf.add_command(band_gap)
minushalf.add_command(create_input)
minushalf.add_command(correct_potfile)
minushalf.add_command(execute)
