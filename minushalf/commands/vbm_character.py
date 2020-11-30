"""
Aims to show how the last conduction band and the
first valence band are composed by the orbitals of each atom.
"""
import click
import pandas as pd
from pathlib import Path
from minushalf.softwares import Vasp
from minushalf.utils import welcome_message, end_message, Orbital, OrbitalType


@click.command()
@click.option('-s',
              '--software',
              type=click.Choice(['VASP'], case_sensitive=False),
              default="VASP",
              show_default=True,
              help="""Specifies the software used to define the
              structure of the file containing the atoms potential.""")
@click.option('-p',
              '--procar-path',
              type=click.Path(exists=True),
              default="PROCAR",
              nargs=1,
              show_default=True,
              help="""Path to PROCAR file. This is only used if the software
              choosed were VASP.""")
@click.option('-e',
              '--eigenval-path',
              default="EIGENVAL",
              type=click.Path(exists=True),
              nargs=1,
              show_default=True,
              help="""Path to EIGENVAL file. This is only used if the software
              choosed were VASP.""")
@click.option(
    '-v',
    '--vasprun-path',
    type=click.Path(exists=True),
    default="vasprun.xml",
    nargs=1,
    show_default=True,
    help="""Path to vasprun.xml file. This is only used if the software
              choosed were VASP.""")
def vbm_character(software: str, procar_path: str, eigenval_path: str,
                  vasprun_path: str) -> None:
    """Uses softwares output files about projections in bands to
    calculate its character. It has to receive a path to an specific file, the list
    of the default names for each software is find bellow:

    VASP: PROCAR, EIGENVAL, vasprun.xml
    """

    welcome_message("minushalf")

    optional_params = {}

    if software == "VASP":
        optional_params = {
            "procar_path": Path(procar_path).__str__(),
            "eigenval_path": Path(eigenval_path).__str__(),
            "vasprun_path": Path(vasprun_path).__str__(),
        }

    softwares = {"VASP": lambda params: Vasp().band_structure(**params)}

    band_structure = softwares[software.upper()](optional_params)
    vbm_projection = band_structure.vbm_projection()
    orbitals = [orbital.__str__() for orbital in Orbital]
    projection_df = (pd.DataFrame.from_dict(vbm_projection,
                                            orient="index",
                                            columns=orbitals))
    normalized_df = projection_df.div(projection_df.sum().sum())
    normalized_df = normalized_df.mul(100).round()

    orbital_type = [str(elem) for elem in OrbitalType]

    def join_orbitals(cols):
        return [
            cat for col in cols for cat in orbital_type if col.startswith(cat)
        ]

    normalized_df = normalized_df.groupby(join_orbitals(projection_df.columns),
                                          axis=1).sum()
    click.echo(normalized_df.to_markdown())

    end_message()
