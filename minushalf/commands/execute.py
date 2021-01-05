"""
Execute command
"""
import os
import sys
import shutil
import click
from loguru import logger
from minushalf.utils import (
    MinushalfYaml,
    BandStructure,
    projection_to_df,
    welcome_message,
    end_message,
)
from minushalf.softwares import (VaspFactory)
from minushalf.corrections import (VaspCorrectionFactory)
from minushalf.data import (Softwares, CorrectionCode)
from minushalf.interfaces import (SoftwaresAbstractFactory)


def get_vbm_projection(factory: SoftwaresAbstractFactory):
    """
    Returns vbm projection
    """
    eigenvalues = factory.get_eigenvalues()
    fermi_energy = factory.get_fermi_energy()
    atoms_map = factory.get_atoms_map()
    num_bands = factory.get_number_of_bands()
    band_projection_file = factory.get_band_projection_class()

    band_structure = BandStructure(eigenvalues, fermi_energy, atoms_map,
                                   num_bands, band_projection_file)
    vbm_projection = band_structure.vbm_projection()
    normalized_df = projection_to_df(vbm_projection)
    return normalized_df


def get_cbm_projection(factory: SoftwaresAbstractFactory):
    """
    Returns cbm projection
    """
    eigenvalues = factory.get_eigenvalues()
    fermi_energy = factory.get_fermi_energy()
    atoms_map = factory.get_atoms_map()
    num_bands = factory.get_number_of_bands()
    band_projection_file = factory.get_band_projection_class()

    band_structure = BandStructure(eigenvalues, fermi_energy, atoms_map,
                                   num_bands, band_projection_file)
    cbm_projection = band_structure.cbm_projection()
    normalized_df = projection_to_df(cbm_projection)
    return normalized_df


def get_atoms_list(factory: SoftwaresAbstractFactory) -> list:
    """
    Returns atoms_list
    """
    atoms_map = factory.get_atoms_map()
    return list(atoms_map.values())


def make_output_file(
    gap: float,
    valence_cuts: dict,
    conduction_cuts: dict = None,
    name: str = "minushalf_results.dat",
) -> None:
    """
    Make output file for execute command
    """
    logger.info("Writing output file")
    with open(name, "w") as file:
        ## Write valence cuts
        file.write("Valence correction cuts:\n")
        for key, value in valence_cuts.items():
            symbol, orbital = key
            cut = value
            file.write("\t({},{}):{:.2f}A\n".format(symbol, orbital, cut))
        file.write(
            "----------------------------------------------------------------")
        ## Write conduction cuts
        if conduction_cuts:
            file.write("Conduction correction cuts:\n")
            for key, value in conduction_cuts.items():
                symbol, orbital = key
                cut = value
                file.write("\t({},{}):{:.2f}A\n".format(symbol, orbital, cut))
        file.write(
            "----------------------------------------------------------------")
        ## Write gap
        file.write("GAP: {:.3}eV".format(gap))


@click.command()
@click.option('--quiet', default=False, is_flag=True)
def execute(quiet: bool):
    """
        Execute command
            Requires:
                minushalf.yaml file : Parameters file
            Return:
                minushalf_results.dat : Results file
    """
    welcome_message("minushalf")

    if quiet:
        logger.remove()
        logger.add(sys.stdout, level="ERROR")
    ## Read yaml file
    logger.info("Reading minushalf.yaml file")
    minushalf_yaml = MinushalfYaml.from_file()
    correction_factory_chooser = {
        Softwares.vasp.value: VaspCorrectionFactory()
    }
    software_factory_chooser = {Softwares.vasp.value: VaspFactory()}

    correction_factory = correction_factory_chooser[minushalf_yaml.software]
    software_factory = software_factory_chooser[minushalf_yaml.software]

    ## Makes abinition calculation
    logger.info("Running ab initio calculations")
    runner = software_factory.get_runner()
    runner.run()

    ## Makes root folder
    logger.info("Make potfiles folder")
    root_folder = "mkpotfiles"
    if os.path.exists(root_folder):
        shutil.rmtree(root_folder)
    os.mkdir(root_folder)
    ## make valence folder
    valence_path = os.path.join(root_folder, "valence")
    os.mkdir(valence_path)

    ## make conduction folder
    conduction_path = os.path.join(root_folder, "conduction")
    os.mkdir(conduction_path)

    ## get vbm projection
    logger.info("Get Vbm and CBM projections")
    vbm_projection = get_vbm_projection(software_factory)
    cbm_projection = get_cbm_projection(software_factory)

    ## get atoms list
    logger.info("Get atoms list")
    atoms = get_atoms_list(software_factory)

    valence_options = {
        "root_folder": valence_path,
        "software_factory": software_factory,
        "runner": runner,
        "minushalf_yaml": minushalf_yaml,
        "vbm_projection": vbm_projection,
        "atoms": atoms
    }
    conduction_options = {
        "root_folder": valence_path,
        "software_factory": software_factory,
        "runner": runner,
        "minushalf_yaml": minushalf_yaml,
        "cbm_projection": cbm_projection,
        "atoms": atoms
    }

    logger.info("Doing corrections")
    if minushalf_yaml.correction.correction_code == CorrectionCode.v.name:

        valence_correction = correction_factory.valence(**valence_options)
        valence_cuts, valence_gap = valence_correction.execute()
        make_output_file(valence_cuts=valence_cuts, gap=valence_gap)

    elif minushalf_yaml.correction.correction_code == CorrectionCode.vf.name:

        valence_fractionary_correction = correction_factory.valence_fractionary(
            **valence_options)
        valence_cuts, valence_gap = valence_fractionary_correction.execute()
        make_output_file(valence_cuts=valence_cuts, gap=valence_gap)

    elif minushalf_yaml.correction.correction_code == CorrectionCode.vc.name:
        valence_correction = correction_factory.valence(**valence_options)
        conduction_correction = correction_factory.conduction(
            **conduction_options)
        valence_cuts, _ = valence_correction.execute()
        conduction_cuts, conduction_gap = conduction_correction.execute()
        make_output_file(valence_cuts=valence_cuts,
                         gap=conduction_gap,
                         conduction_cuts=conduction_cuts)

    elif minushalf_yaml.correction.correction_code == CorrectionCode.vfc.name:

        valence_fractionary_correction = correction_factory.valence_fractionary(
            **valence_options)
        conduction_correction = correction_factory.conduction(
            **conduction_options)
        valence_cuts, _ = valence_fractionary_correction.execute()
        conduction_cuts, conduction_gap = conduction_correction.execute()
        make_output_file(valence_cuts=valence_cuts,
                         gap=conduction_gap,
                         conduction_cuts=conduction_cuts)

    elif minushalf_yaml.correction.correction_code == CorrectionCode.vfcf.name:

        valence_fractionary_correction = correction_factory.valence_fractionary(
            **valence_options)
        conduction_fractionary_correction = correction_factory.conduction_fractionary(
            **conduction_options)
        valence_fractionary_correction.execute()
        conduction_fractionary_correction.execute()

    elif minushalf_yaml.correction.correction_code == CorrectionCode.vcf.name:
        valence_correction = correction_factory.valence(**valence_options)
        conduction_fractionary_correction = correction_factory.conduction_fractionary(
            **conduction_options)
        valence_cuts, _ = valence_correction.execute()
        conduction_cuts, conduction_gap = conduction_fractionary_correction.execute(
        )
        make_output_file(valence_cuts=valence_cuts,
                         gap=conduction_gap,
                         conduction_cuts=conduction_cuts)
    end_message()
