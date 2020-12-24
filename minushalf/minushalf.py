"""
Definition of the minushalf CLI
"""
import click
from minushalf.commands import (run_atomic, occupation, vbm_character,
                                cbm_character, band_character, create_input)


@click.group()
def minushalf():
    """
    Group cli minushalf
    """


minushalf.add_command(run_atomic)
minushalf.add_command(occupation)
minushalf.add_command(vbm_character)
minushalf.add_command(cbm_character)
minushalf.add_command(band_character)
minushalf.add_command(create_input)
