"""
Atomic Correct potential file
"""
from __future__ import annotations
import sys
import click
from numpy import arange
from loguru import logger
from minushalf.utils import welcome_message, end_message
from minushalf.corrections import Corrections


@click.command()
@click.option('--quiet', default=False, is_flag=True)
@click.option('-p',
              '--pot_path',
              type=click.Path(exists=True),
              show_default=True,
              help="""Path relative to the potential files.
              If this variable is not provided, the default
              name for each software will be used, as described below:

              VASP: POTCAR""")
@click.option('-s',
              '--software',
              type=click.Choice(['VASP'], case_sensitive=False),
              default="VASP",
              show_default=True,
              help="""Specifies the software used to define the
              structure of the file containing the atoms potential.""")
@click.option('-c',
              '--correction',
              type=click.Choice(['VALENCE', 'CONDUCTION'],
                                case_sensitive=False),
              default="VALENCE",
              show_default=True,
              help="""Indicates whether the correction should be
                      made in the last valence band, the first conduction band or both."""
              )
@click.option('-C',
              '--cut',
              type=str,
              nargs=1,
              help="""
              distance value used to cut the potential generated artificially
              by fractional atomic occupation , it can be passed in two ways:

              unique value : float or integer

              \b
              range:  begin(float|integer):pass(float|integer):end(float|integer)

              """)
@click.option('-a',
              '--amplitude',
              type=click.FloatRange(0.0, 20.0),
              default=1.0,
              show_default=True,
              nargs=1,
              help="""Multiplying factor to be used to
              correct the artificially generated potential.""")
def correct_potfile(
    quiet: bool,
    pot_path: str,
    software: str,
    correction: str,
    cut: str,
    amplitude: float,
) -> None:
    """Generate occupied atomic potential file used for ab initio calculations.

    Requires:

        VTOTAL.ae

        VTOTAL_OCC

        POTENTIAL_FILE(POTCAR)

        INP_OCC


    Generates:

        POTFILEcut${CUT_VALUE} (If amplitude is equal to 1.0)

        POTFILEcut${CUT_VALUE}A${AMP_VALUE} (If amplitude is different to 1.0)

    """
    welcome_message("minushalf")

    if quiet:
        logger.remove()
        logger.add(sys.stdout, level="ERROR")

    if correction == 'CONDUCTION':
        amplitude *= -1

    apply_correction(software, amplitude, pot_path, cut)

    end_message()


def apply_correction(software: str, amplitude: float, potential_file_path: str,
                     cut: str) -> None:
    """
    Apply corrections on potencial file
    """

    logger.info("Applying correction")

    correction_factory = Corrections()
    softwares = {
        "VASP":
        lambda amplitude, pot_file, new_cut: correction_factory.vasp(
            new_cut, amplitude, pot_file)
        if pot_file else correction_factory.vasp(new_cut, amplitude)
    }

    cut_numbers = parse_cut(cut)

    for new_cut in cut_numbers:
        logger.info("Correcting POTFILE for cut = {:.3} ".format(new_cut))
        program = softwares[software.upper()](amplitude, potential_file_path,
                                              new_cut)

        program.correction()


def parse_cut(cut: str) -> list:
    """
    Parse cut in a list of numbers.

    Args:
        cut (str): Cut energy to be used in the program, it can be
                   passed in two ways:

                    unique value : float or integer
                    range:  begin(float|integer):pass(float|integer):end(float|integer)
    Returns:

        cut_numbers (list): Permited values of cut.
    """
    if not cut:
        raise ValueError('A valeu of cut energy must be provided.')

    cuts = cut.split(':')

    try:
        cut_number = [float(element) for element in cuts]
    except:
        raise ValueError("Invalid Input.")

    if len(cut_number) == 1:
        return cut_number
    elif len(cut_number) != 3:
        raise ValueError()

    return arange(cut_number[0], cut_number[2], cut_number[1])
