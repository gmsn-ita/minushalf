"""
Aims to show how the last conduction band and the
first valence band are composed by the orbitals of each atom.
"""
import click
from minushalf.softwares.softwares import Softwares, get_software_factory
from minushalf.utils.cli_messages import welcome_message, end_message
from minushalf.utils.projection_to_df import projection_to_df
from minushalf.utils.band_structure import BandStructure
from minushalf.utils.software_output import get_output_filenames



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
@click.option('-n',
              '--input-name',
              type=click.Path(),
              nargs=1,
              default=None,
              help="""Path to the input file for pw.x or projwfc.x. Required for QE to determine
            <prefix> of output files."""
)
def cbm_character(software: str, base_path: str, input_name: str, indirect: bool) -> None:
    """Uses output files from softwares that perform ab initio calculations to
      discover the first conduction band (CBM) and extract, in percentage, its
      character corresponding to each orbital type (s, p, d, ... ). The
      names of the files required for each software are listed below, it is
      worth mentioning that their names cannot be modified.

    VASP: PROCAR, EIGENVAL, vasprun.xml
    QE:   {prefix}.xml, {prefix}.save/, {prefix}.pdos_atm* (requires --input-name)
    """

    welcome_message("minushalf")

    filenames = get_output_filenames(software.upper(), input_name, base_path)
    factory = get_software_factory(software.upper())

    eigenvalues          = factory.get_eigenvalues(
                               filename=filenames["eigenvalues"],
                               base_path=base_path)
    fermi_energy         = factory.get_fermi_energy(
                               filename=filenames["fermi_energy"],
                               base_path=base_path)
    atoms_map            = factory.get_atoms_map(
                               filename=filenames["atoms_map"],
                               base_path=base_path)
    num_bands            = factory.get_number_of_bands(
                               filename=filenames["number_of_bands"],
                               base_path=base_path)
    band_projection_file = factory.get_band_projection_class(
                               filename=filenames["band_projection"],
                               base_path=base_path)

    band_structure = BandStructure(eigenvalues, fermi_energy, atoms_map,
                                   num_bands, band_projection_file)
    cbm_index = band_structure.cbm_index(is_indirect=indirect)
    click.echo(f"CBM: Kpoint {cbm_index[0]}, band {cbm_index[1]}")
    cbm_projection = band_structure.cbm_projection(is_indirect=indirect)
    normalized_df = projection_to_df(cbm_projection)
    click.echo(normalized_df.to_markdown())

    end_message()
