"""
Aims to show how the last conduction band and the
first valence band are composed by the orbitals of each atom.
"""
import click
from minushalf.data import Softwares
from minushalf.softwares import VaspFactory
from minushalf.utils import (welcome_message, end_message, projection_to_df,
                             BandStructure)


@click.command()
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
def cbm_character(software: str, base_path: str) -> None:
    """Uses output files from softwares that perform ab initio calculations to
      discover the first conduction band (CBM) and extract, in percentage, its
      character corresponding to each orbital type (s, p, d, ... ). The
      names of the files required for each software are listed below, it is
      worth mentioning that their names cannot be modified.

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
    cbm_projection = band_structure.cbm_projection()
    normalized_df = projection_to_df(cbm_projection)
    click.echo(normalized_df.to_markdown())

    end_message()
