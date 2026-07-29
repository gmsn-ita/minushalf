"""
Makes fractional occupation on INP file
"""
import sys
import click
from loguru import logger
from minushalf.softwares.softwares import Softwares
from minushalf.utils.calculation_code import CalculationCode
from minushalf.utils.exchange_correlation import ExchangeCorrelation
from minushalf.io.input_file import (InputFile)
from minushalf.utils.cli_messages import (welcome_message, end_message)


@click.command()
@click.argument('chemical_symbol', type=str, nargs=1)
@click.option(
    '-s',
    '--software',
    type=click.Choice(Softwares.to_list(), case_sensitive=False),
    default=Softwares.vasp.value,
    show_default=True,
    help="""Specifies the software used to perform ab initio calculations.""")
@click.option('-e',
              '--exchange_correlation_code',
              type=click.Choice(ExchangeCorrelation.to_list(),
                                case_sensitive=False),
              nargs=1,
              default=ExchangeCorrelation.pb.value,
              show_default=True,
              help="""
            Represents the functional of exchange and correlation, it can assume the following values:


              ca: Ceperley-Alder (VASP only)

              wi: Wigner

              hl: Hedin-Lundqvist

              gl: Gunnarson-Lundqvist 

              bh: Von Barth-Hedin (VASP only)

              pb: PBE scheme by Perdew, Burke, and Ernzerhof

              rp: RPBE scheme by Hammer, Hansen, and Norskov

              rv: revPBE scheme by Zhang and Yang

              bl: BLYP (Becke-Lee-Yang-Parr) scheme

              pz: Perdew-Zunger LDA (Quantum ESPRESSO only)

              """)
# This is the "dft" tag in Quantum ESPRESSO
@click.option('-c',
              '--calculation_code',
              type=click.Choice(CalculationCode.to_list(),
                                case_sensitive=False),
              nargs=1,
              default=CalculationCode.ae.value,
              show_default=True,
              help="""Represents calculation code,it can
              assume the following values:

              ae: All electrons
              """)
@click.option(
    '-m',
    "--maximum_iterations",
    type=click.IntRange(0, 10000),
    nargs=1,
    default=100,
    show_default=True,
    help="""Maximum number of iterations performed by the atomic program""")
@click.option('-f',
              "--filename",
              type=str,
              nargs=1,
              default='INP',
              show_default=True,
              help="""Name of the created file""")
@click.option('-r',
              "--cut",
              type=float,
              nargs=1,
              default='0.0',
              show_default=True,
              help="""CUT parameter""")
@click.option('--quiet', default=False, is_flag=True)
def create_input(software: str, cut: float, chemical_symbol: str, exchange_correlation_code: str,
                 calculation_code: str, maximum_iterations: int, filename: str,
                 quiet: bool):
    """
    Create the input file for the run-atomic command.


    Requires:

        CHEMICAL_SYMBOL: Chemical symbol of the atom (H, He, Na, Li...). Check the list
                         of available atoms in the docs


    Returns:

        INP: The input file for run-atomic command
    """

    welcome_message("minushalf")

    if quiet:
        logger.remove()
        logger.add(sys.stdout, level="ERROR")

    input_file = InputFile.minimum_setup(chemical_symbol.capitalize(),
                                         exchange_correlation_code,
                                         maximum_iterations, calculation_code,
                                         software.upper(), cut)
    logger.info("Creating INP file")  

    input_file.to_file(filename)

    end_message()
