"""
Definition of the minushalf CLI
"""
import click
from minushalf.commands import (run_atomic, occupation, correct_potfile,
                                vbm_character, cbm_character, band_character)


@click.group()
def minushalf():
    """
    Group cli minushalf
    """


minushalf.add_command(run_atomic)
minushalf.add_command(occupation)
minushalf.add_command(correct_potfile)
minushalf.add_command(vbm_character)
minushalf.add_command(cbm_character)
minushalf.add_command(band_character)
