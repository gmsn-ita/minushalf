"""
Aims to show how the last valence band are composed by the orbitals of each atom.
"""
import click
from minushalf.softwares.softwares import Softwares, get_software_factory
from minushalf.utils.cli_messages import welcome_message, end_message
from minushalf.utils.projection_to_df import projection_to_df
from minushalf.utils.band_structure import BandStructure


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
@click.option('-i',
              '--indirect',
              type=bool,
              nargs=1,
              is_flag=True,
              help="""Calculate indirect band gap.""")
def vbm_character(software: str, base_path: str, indirect: bool) -> None:
    """Uses output files from softwares that perform ab initio calculations to discover the
       last valence band (VBM) and extract, in percentage, its character corresponding to each
       orbital type (s, p, d, ... ). The names of the files required for each
       software are listed below, it is worth mentioning that their names cannot be modified.

    VASP: PROCAR, EIGENVAL, vasprun.xml
    """

    welcome_message("minushalf")

    factory = get_software_factory(software.upper())

    eigenvalues = factory.get_eigenvalues(base_path=base_path)
    fermi_energy = factory.get_fermi_energy(base_path=base_path)
    atoms_map = factory.get_atoms_map(base_path=base_path)
    num_bands = factory.get_number_of_bands(base_path=base_path)
    band_projection_file = factory.get_band_projection_class(
        base_path=base_path)

    band_structure = BandStructure(eigenvalues, fermi_energy, atoms_map,
                                   num_bands, band_projection_file)
    vbm_index = band_structure.vbm_index(is_indirect=indirect)
    click.echo(f"VBM: Kpoint {vbm_index[0]}, band {vbm_index[1]}")
    vbm_projection = band_structure.vbm_projection(is_indirect=indirect)
    normalized_df = projection_to_df(vbm_projection)

    click.echo(normalized_df.to_markdown())

    end_message()
