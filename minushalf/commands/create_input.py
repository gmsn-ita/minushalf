"""
Makes fractional occupation on INP file
"""
import sys
import click
from loguru import logger
from minushalf.data import (CalculationCode, ExchangeCorreltion)
from minushalf.utils import (welcome_message, end_message, InputFile)


@click.command()
@click.argument('chemical_symbol', type=str, nargs=1)
@click.option('-e',
              '--exchange_correlation_code',
              type=click.Choice(ExchangeCorreltion.to_list(),
                                case_sensitive=False),
              nargs=1,
              default=ExchangeCorreltion.pb.value,
              show_default=True,
              help="""
            Represents the functional of exchange and correlation,it can assume the following values:


              ca: Ceperley-Alder

              wi: Wigner

              hl: Hedin-Lundqvist

              gl: Gunnarson-Lundqvist

              bh: Von Barth-Hedin

              pb: PBE scheme by Perdew, Burke, and Ernzerhof

              rp: RPBE scheme by Hammer, Hansen, and Norskov

              rv: revPBE scheme by Zhang and Yang

              bl: BLYP (Becke-Lee-Yang-Parr) scheme

              """)
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
@click.option('--quiet', default=False, is_flag=True)
def create_input(chemical_symbol: str, exchange_correlation_code: str,
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
                                         maximum_iterations, calculation_code)
    logger.info("Creating INP file")

    input_file.to_file(filename)

    end_message()
