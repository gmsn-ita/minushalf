"""
Aims to show how the last conduction band and the
first valence band are composed by the orbitals of each atom.
"""
import click
from minushalf.projected_wave_function import ProjectedWaveFunction
from minushalf.utils import welcome_message, end_message


@click.command()
@click.option(
    "-v",
    "--vbm",
    "report",
    flag_value="valence",
    help="Print only vbm report.",
)
@click.option(
    "-c",
    "--cbm",
    "report",
    flag_value="conduction",
    help="Print only cbm report.",
)
@click.option(
    "-b",
    "--both",
    "report",
    flag_value="both",
    help="Print vbm and cbm reports.",
    default=True,
    show_default=True
)
@click.option('-s', '--software',
              type=click.Choice(['VASP'],
                                case_sensitive=False),
              default="VASP",
              show_default=True,
              help="""Specifies the software used to define the
              structure of the file containing the atoms potential.""")
@click.argument("path", type=click.Path(exists=True), nargs=1)
def band_character(path: click.Path, report: str, software: str) -> None:
    """Uses softwares output files about projections in bands to
    calculate its character. It has to receive a path to an specific file, the list
    of the default names for each software is find bellow:

    VASP: vasprun.xml
    """

    welcome_message("minushalf")
    factory = ProjectedWaveFunction()

    softwares = {
        "VASP": lambda file_path: factory.vasp(file_path) if file_path else factory.vasp(),
    }

    soft_projected_wave_function = softwares[software](path.__str__())

    if report == "valence" or report == "both":
        click.echo(click.style("Projection VBM:\n", bold=True))
        band_result_vbm = soft_projected_wave_function.character_vbm()
        click.echo(click.style(band_result_vbm.to_markdown(), bold=True))

    if report == "conduction" or report == "both":
        click.echo(click.style("Projection CBM:\n", bold=True))
        band_result_cbm = soft_projected_wave_function.character_cbm()
        click.echo(click.style(band_result_cbm.to_markdown(), bold=True))

    end_message()
