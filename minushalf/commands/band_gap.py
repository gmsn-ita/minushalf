"""
Command to read band gap
"""

import click
from minushalf.softwares import VaspFactory
from minushalf.data import Softwares
from minushalf.utils import (welcome_message, end_message, BandStructure)


@click.command()
@click.option('-s',
              '--software',
              type=click.Choice(Softwares.to_list(), case_sensitive=False),
              default=Softwares.vasp.value,
              show_default=True,
              help="""Specifies the software used to define the
              structure of the file containing the atoms potential.""")
@click.option(
    '-b',
    '--base-path',
    type=click.Path(),
    nargs=1,
    help="""folder where the necessary files for the calculation are located"""
)
def band_gap(software: str, base_path: str) -> None:
    """Uses softwares output files about projections in bands to
    calculate its character. It has to receive a path to an specific file, the list
    of the default names for each software is find bellow:

    VASP: PROCAR, EIGENVAL, vasprun.xml
    """

    welcome_message("minushalf")

    softwares = {"VASP": VaspFactory()}

    factory = softwares[software.upper()]

    eigenvalues = factory.get_eigenvalues(base_path=base_path)
    fermi_energy = factory.get_fermi_energy(base_path=base_path)
    atoms_map = factory.get_atoms_map(base_path=base_path)
    num_bands = factory.get_number_of_bands(base_path=base_path)
    band_projection_file = factory.get_band_projection_class(
        base_path=base_path)

    band_structure = BandStructure(eigenvalues, fermi_energy, atoms_map,
                                   num_bands, band_projection_file)

    gap_report = band_structure.band_gap()
    click.echo(gap_report["vbm"])
    click.echo(gap_report["cbm"])
    click.echo("Gap: {:.3f}eV".format(gap_report["gap"]))
    end_message()
