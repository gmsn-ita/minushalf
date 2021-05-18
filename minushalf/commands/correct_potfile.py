"""
Atomic Correct potential file
"""
import sys
import click
from loguru import logger
from minushalf.softwares import VaspFactory
from minushalf.data import Softwares
from minushalf.utils import (welcome_message, end_message, Vtotal,
                             AtomicPotential, parse_cut)


@click.command()
@click.option('--quiet', default=False, is_flag=True)
@click.option('-b',
              '--base_potfile_path',
              type=click.Path(),
              help="""Path to the folder containing the potential file""")
@click.option('-v',
              '--vtotal_path',
              type=click.Path(exists=True),
              default="VTOTAL.ae",
              show_default=True,
              help="""Path relative to the potential
              file generated by the atomic program for the
              atom with all electrons.""")
@click.option('-o',
              '--vtotal_occupied_path',
              type=click.Path(exists=True),
              default="VTOTAL_OCC",
              show_default=True,
              help="""Path relative to the potential
              file generated by the atomic program for the
              occupied atom.""")
@click.option(
    '-s',
    '--software',
    type=click.Choice(Softwares.to_list(), case_sensitive=False),
    default=Softwares.vasp.value,
    show_default=True,
    help="""Specifies the software used to make ab initio calculations.""")
@click.option('-c',
              '--correction',
              type=click.Choice(['VALENCE', 'CONDUCTION'],
                                case_sensitive=False),
              default="VALENCE",
              show_default=True,
              help="""Indicates whether the correction should be
                      made in the valence band or the conduction band.""")
@click.option('-C',
              '--cut',
              type=str,
              nargs=1,
              default="2.0",
              show_default=True,
              help="""
              distance value used to cut the potential generated artificially
              by fractional atomic occupation , it can be passed in two ways:

              unique value : float or integer. Ex: 1.0

              \b
              range:  begin(float|integer):pass(float|integer):end(float|integer). Ex: 1.0:0.1:2.0

              """)
@click.option(
    '-a',
    '--amplitude',
    type=float,
    default=1.0,
    show_default=True,
    nargs=1,
    help=
    """Scaling factor to be used to correct the artificially generated potential.
               In the vast majority of cases, the amplitude value is 1.0. However, there are some
               special cases where this value needs to be adjusted. Therefore, we recommend that
                you do not change this value unless you know exactly what you are doing"""
)
def correct_potfile(
    quiet: bool,
    base_potfile_path: str,
    vtotal_path: str,
    vtotal_occupied_path: str,
    software: str,
    correction: str,
    cut: str,
    amplitude: float,
) -> None:
    """Generate the occupied atomic potential file used for ab initio calculations.

    Requires:

        VTOTAL.ae: potential of the atom with all electrons

        VTOTAL_OCC: potential of the occupied atom

        INP_OCC: Input file for the run-atomic command of the occupied atom

        The command also needs the potential files used by the chosen software:

            VASP: POTCAR (This name can't be changed)


    Generates:

        POTFILEcut${CUT_VALUE} (If amplitude is equal to 1.0)

        POTFILEcut${CUT_VALUE}A${AMPLITUDE_VALUE} (If amplitude is different from 1.0)

    """
    welcome_message("minushalf")

    if quiet:
        logger.remove()
        logger.add(sys.stdout, level="ERROR")

    softwares = {"VASP": VaspFactory()}
    factory = softwares[software.upper()]

    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occupied_path)
    potential_file = factory.get_potential_class(base_path=base_potfile_path)

    atomic_potential = AtomicPotential(vtotal, vtotal_occ, potential_file)

    cut_numbers = parse_cut(cut)

    is_conduction = False
    if correction.upper() == "CONDUCTION":
        is_conduction = True

    for new_cut in cut_numbers:
        logger.info("Correcting POTFILE for cut = {:.3} ".format(new_cut))
        new_potential = atomic_potential.correct_potential(
            new_cut, amplitude, is_conduction)
        atomic_potential.correct_file(new_potential, new_cut, amplitude,
                                      is_conduction)

    end_message()
