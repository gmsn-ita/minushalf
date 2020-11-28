"""
Definition of the minushalf CLI
"""
import click
from .commands import (run_atomic, band_character, occupation, correct_potfile)


@click.group()
def minushalf():
    """
    Group cli minushalf
    """


minushalf.add_command(run_atomic)
minushalf.add_command(band_character)
minushalf.add_command(occupation)
minushalf.add_command(correct_potfile)
