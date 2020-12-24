"""
Run atomic program
"""
import os
import sys
import click
from loguru import logger
from minushalf.utils import welcome_message, end_message
from minushalf import atomic_program


@click.command()
@click.option('--quiet', default=False, is_flag=True)
def run_atomic(quiet: bool):
    """
    Run atomic program.

    Requires:

        INP: The input file for the calculation.

    Returns:

        INP.ae: A copy of the input file for the calculation.


        VTOTAL.ae: Contains the atom potential.


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

    logger.info("Run atomic program")

    try:
        atomic_program.run()
    except:
        raise Exception('Problems in atomic program execution')

    logger.info("Atomic program finished execution.")

    if not os.path.exists('./VTOTAL1'):
        raise FileNotFoundError("VTOTAL0 not found")

    logger.info("Changing VTOTAL1 to VTOTAL.ae")
    os.rename("VTOTAL1", "VTOTAL.ae")

    end_message()
