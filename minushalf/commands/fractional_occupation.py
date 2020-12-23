"""
Makes fractional occupation on INP file
"""
import os
import sys
import click
import numpy as np
from loguru import logger
from minushalf.utils import (welcome_message, end_message, InputFile)
from minushalf import atomic_program


@click.command()
@click.argument('orbital_quantum_number', type=str, nargs=1)
@click.argument('occupation_percentual', type=str, nargs=1, default="100")
@click.option('--quiet', default=False, is_flag=True)
def occupation(orbital_quantum_number: str, occupation_percentual: str,
               quiet: bool):
    """
    Perform fractional occupation on the atom and generate the pseudopotential for this occupation.
    The occupation can subtract any fraction of the electron between 0 and 0.5, half occupation is the default.

    Requires:

        ORBITAL_QUANTUM_NUMBER: A string that defines the orbital(s) in which the occupation will be made,
        it can assume four values: (0: s | 1: p | 2: d | 3: f). if going to pass multiple orbitals,
        pass a string with numbers separated by commas : ("0,1,2,3")


        OCCUPATION_PERCENTUAL: A string that defines percentual of half eletron to be used in the occupation.
        The default is 100%, wich states for 0.5e. For multiple occupations in different orbitals, pass a string
        separated by commas ("100,50,40,100"). For simplicity, to avoid the excessive repetition of the number
        100, just replace the number with * ("*,30,*"). If this argument is not used, the occupation of 
        half electron will be made for all orbitals


        INP: A copy of the input file used in ATOM program

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

    input_file = InputFile.from_file()
    logger.info("Adding minus one half electron correction on INP")

    try:
        quantum_numbers = np.array(orbital_quantum_number.split(","),
                                   dtype=np.int)
    except ValueError as wrong_input:
        raise ValueError(
            "Invalid value for secondary quantum number") from wrong_input

    for quantum_number in quantum_numbers:
        if quantum_number < 0 or quantum_number > 3:
            raise ValueError("Invalid value for secondary quantum number")

    percentuals = np.ones(len(quantum_numbers)) * 100

    for index, percentual in enumerate(occupation_percentual.split(",")):
        if percentual != "*":
            if float(percentual) < 0 or float(percentual) > 100:
                raise ValueError("Invalid value for occupation percentual")
            percentuals[index] = float(percentual)

    for quantum_number, percentual in zip(quantum_numbers, percentuals):
        input_file.electron_occupation(0.5 * (percentual / 100),
                                       quantum_number)

    os.rename('INP', 'INP.ae')
    input_file.to_file()

    logger.info("Run atomic program")
    try:
        atomic_program.run()
    except Exception as program_fail:
        raise Exception('Problems in atomic program') from program_fail

    logger.info("Atomic program finished execution.")

    if not os.path.exists('./VTOTAL1'):
        raise FileNotFoundError("Problems in INP file generation")
    logger.info("Changing VTOTAL1 to VTOTAL_OCC")
    os.rename("VTOTAL1", "VTOTAL_OCC")

    logger.info("Changing INP to INP_OCC")
    os.rename('INP', 'INP_OCC')

    end_message()
