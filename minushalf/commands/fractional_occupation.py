"""
Makes fractional occupation on INP file
"""
from __future__ import annotations
import os
import sys
import click
from loguru import logger
from minushalf.utils import welcome_message, end_message
from minushalf.atomic import InputFile
from minushalf import atomic_program


@click.command()
@click.argument('occupation_percentual', type=click.FloatRange(0, 100),
                nargs=1,
                default=100)
@click.option('--quiet', default=False, is_flag=True)
def occupation(occupation_percentual: float, quiet: bool):
    """
    Perform fractional occupation on the atom and generate the potential for this occupation.
    The occupation can contain any fraction of the electron between 0 and 0.5, half occupation is the default.

    Requires:


        OCCUPATION_PERCENTUAL: The percentual of half eletron to be used in the occupation.
        The default is 100%, wich states for 0.5e.


        INP: A copy of the input file for the calculation.

    Returns:

        INP_OCC : Input file modified for fractional occupation


        INP.ae: A copy of the input file for the calculation.


        VTOTAL_OCC: Contains the atom potential for fractional occupation.


        OUT: Contains detailed information about the run.


        AECHARGE: Contains in four columns values of r, the “up” and “down” parts of the total
        charge density, and the total core charge density (the charges multiplied by 4πr 2 ).


        CHARGE: is exactly identical to AECHARGE and is generated for backwards compatibility.


        RHO: Like CHARGE, but without the 4πr 2 factor


        AEWFNR0...AEWFNR3: All-electron valence wavefunctions as function of radius, for s, p, d,
        and f valence orbitals (0, 1, 2, 3, respectively — some channels might not be available).
        They include a factor of r, the s orbitals also going to zero at the origin.
    """

    welcome_message("minushalf")

    if quiet:
        logger.remove()
        logger.add(sys.stdout, level="ERROR")

    input_file = InputFile()
    logger.info("Adding minus one half electron correction on INP")
    input_file.electron_occupation(0.5*(occupation_percentual/100))

    os.rename('INP', 'INP.ae')
    input_file.to_file()

    logger.info("Run atomic program")
    try:
        atomic_program.run()
    except:
        raise Exception('Problems in atomic program')

    logger.info("Atomic program finished execution.")

    if not os.path.exists('./VTOTAL1'):
        raise FileNotFoundError("Problems in INP file generation")
    logger.info("Changing VTOTAL1 to VTOTAL_OCC")
    os.rename("VTOTAL1", "VTOTAL_OCC")

    logger.info("Changing INP to INP_OCC")
    os.rename('INP', 'INP_OCC')

    logger.info("Apply correction on potential file")

    end_message()
