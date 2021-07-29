"""
Aims to show how the band in a specific k-point is
composed by the orbitals of each atom.
"""
import click
from minushalf.data import Softwares
from minushalf.softwares import Vasp
from minushalf.utils import (welcome_message, end_message, projection_to_df,
                             BandStructure)


@click.command()
@click.argument('kpoint', nargs=1, type=int)
@click.argument('band', nargs=1, type=int)
@click.option(
    '-s',
    '--software',
    type=click.Choice(Softwares.to_list(), case_sensitive=False),
    default=Softwares.vasp.value,
    show_default=True,
    help="""Specifies the software used to perform ab initio calculations.""")
@click.option('-b',
              '--base-path',
              type=click.Path(),
              nargs=1,
              help="""Path to folder where the relevant files are located.""")
def band_character(kpoint: int, band: int, software: str,
                   base_path: str) -> None:
    """Uses output files from softwares that perform ab initio calculations to
      read projections in a specific kpoint band and extract, in percentage, its
      character corresponding to each orbital type (s, p, d, ... ). The
      names of the files required for each software are listed below, it is
      worth mentioning that their names cannot be modified.

    VASP: PROCAR, EIGENVAL, vasprun.xml
    """

    welcome_message("minushalf")

    softwares = {"VASP": Vasp()}

    factory = softwares[software.upper()]

    eigenvalues = factory.get_eigenvalues(base_path=base_path)
    fermi_energy = factory.get_fermi_energy(base_path=base_path)
    atoms_map = factory.get_atoms_map(base_path=base_path)
    num_bands = factory.get_number_of_bands(base_path=base_path)
    band_projection_file = factory.get_band_projection_class(
        base_path=base_path)

    band_structure = BandStructure(eigenvalues, fermi_energy, atoms_map,
                                   num_bands, band_projection_file)
    band_projection = band_structure.band_projection(kpoint, band)
    normalized_df = projection_to_df(band_projection)

    click.echo(normalized_df.to_markdown())

    end_message()
