"""
Execute command
"""
import os
import sys
import shutil
from collections import OrderedDict
import numpy as np
import pandas as pd
import click
from loguru import logger
from minushalf.utils import (
    MinushalfYaml,
    BandStructure,
    projection_to_df,
    welcome_message,
    end_message,
    make_minushalf_results,
    get_fractionary_corrections_indexes,
    get_simple_corrections_indexes,
)
from minushalf.softwares import (VaspFactory)
from minushalf.corrections import (VaspCorrection)
from minushalf.data import (Softwares, CorrectionCode, CorrectionDefaultParams)
from minushalf.interfaces import (SoftwaresAbstractFactory)


def get_vbm_projection(factory: SoftwaresAbstractFactory) -> pd.DataFrame:
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


def get_cbm_projection(factory: SoftwaresAbstractFactory) -> pd.DataFrame:
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
    atoms = [atoms_map[key] for key in sorted(atoms_map)]
    return list(OrderedDict.fromkeys(atoms))


def overwrite_band_projection(new_values: list,
                              band_projection: pd.DataFrame) -> pd.DataFrame:
    """
    Overwrite values in VBM or CBM band projection

        Args:
            new_values (list): Arguments passed in overwrite_vbm or overwrite_cbm.
            band_projection (pd.Dataframe): Dataframe with the value of projections
                                            of VBM or CBM.
        Returns:
            modified_band_projection (pd.Dataframe): Band projection data frame
                                                    with the values overwrited.
    """
    for case in new_values:
        atom = case[0].capitalize()
        orbital = case[1].lower()
        projection = int(case[2])
        band_projection[orbital][atom] = projection
    return band_projection


@click.command()
@click.option('--quiet', default=False, is_flag=True)
def execute(quiet: bool):
    """

        Uses the Nelder-Mead method to find
        the optimal values for the CUT(S) and,
        finally, find the corrected Gap value.
        This command uses external software to
        perform ab initio calculations, so it must
        be installed in order to perform the command.
        Check the docs for an list of the softwares supported
        by the CLI.


            Requires:


                minushalf.yaml : Parameters file. Check the docs
                                 for a more detailed description.

                ab_initio_files: Files needed to perform the ab initio calculations.
                                 They must be in the same directory as the input
                                 file minushalf.yaml

                potential_folder: Folder with the potential files for each atom in
                                  the crystal. The files must be named in the following pattern
                                  ${POTENTIAL_FILE_NAME}.${LOWERCASE_CHEMICAL_SYMBOL}

            Returns:

                minushalf_results.dat : File that contains the optimal
                                        values of the cuts and the final
                                        value of the Gap.

                corrected_valence_potfiles: Potential files corrected with opti-mum valence cuts.

                corrected_conduction_potfiles: Potential files corrected with optimum conduction cuts.
    """
    welcome_message("minushalf")

    if quiet:
        logger.remove()
        logger.add(sys.stdout, level="ERROR")
    ## Read yaml file
    logger.info("Reading minushalf.yaml file")
    minushalf_yaml = MinushalfYaml.from_file()
    correction_factory_chooser = {Softwares.vasp.value: VaspCorrection}
    software_factory_chooser = {Softwares.vasp.value: VaspFactory()}

    correction = correction_factory_chooser[minushalf_yaml.software]
    software_factory = software_factory_chooser[minushalf_yaml.software]

    ## Makes abinition calculation
    logger.info("Running ab initio calculations")
    software_configurations = minushalf_yaml.software_configurations
    runner = software_factory.get_runner(**software_configurations)
    runner.run()

    ## Makes root folder
    logger.info("Make potfiles folder")
    root_folder = "mkpotfiles"
    if os.path.exists(root_folder):
        shutil.rmtree(root_folder)
    os.mkdir(root_folder)

    ## get vbm projection
    logger.info("Get Vbm and CBM projections")
    vbm_projection = get_vbm_projection(software_factory)
    cbm_projection = get_cbm_projection(software_factory)

    ### Overwrite band projections
    if len(minushalf_yaml.correction[
            CorrectionDefaultParams.overwrite_vbm.name]) > 0:
        logger.warning(
            "You're changing directly the band character. This is not recommendend unless you know exactly what are you doing."
        )
        vbm_projection = overwrite_band_projection(
            minushalf_yaml.correction[
                CorrectionDefaultParams.overwrite_vbm.name], vbm_projection)
    if len(minushalf_yaml.correction[
            CorrectionDefaultParams.overwrite_cbm.name]) > 0:
        logger.warning(
            "You're changing directly the band character. This is not recommendend unless you know exactly what are you doing."
        )
        cbm_projection = overwrite_band_projection(
            minushalf_yaml.correction[
                CorrectionDefaultParams.overwrite_cbm.name], cbm_projection)

    ## get atoms list
    logger.info("Get atoms list")
    atoms = get_atoms_list(software_factory)

    ## amplitude logger
    if not np.isclose(
            minushalf_yaml.correction[CorrectionDefaultParams.amplitude.name],
            CorrectionDefaultParams.amplitude.value):
        logger.warning(
            "Amplitude value is different from 1.0. This is not recommended unless you know exactly what you are doing."
        )

    valence_options = {
        "root_folder": root_folder,
        "software_factory": software_factory,
        "runner": runner,
        "minushalf_yaml": minushalf_yaml,
        "band_projection": vbm_projection,
        "atoms": atoms,
        "is_conduction": False
    }
    conduction_options = {
        "root_folder": root_folder,
        "software_factory": software_factory,
        "runner": runner,
        "minushalf_yaml": minushalf_yaml,
        "band_projection": cbm_projection,
        "atoms": atoms,
        "is_conduction": True
    }

    logger.info("Doing corrections")
    if minushalf_yaml.correction["correction_code"] == CorrectionCode.v.name:
        valence_options["correction_indexes"] = get_simple_corrections_indexes(
            vbm_projection)
        valence_correction = correction(**valence_options)
        valence_cuts, valence_gap = valence_correction.execute()
        make_minushalf_results(valence_cuts=valence_cuts, gap=valence_gap)

    elif minushalf_yaml.correction[
            "correction_code"] == CorrectionCode.vf.name:

        valence_options[
            "correction_indexes"] = get_fractionary_corrections_indexes(
                vbm_projection,
                treshold=minushalf_yaml.
                correction["fractionary_valence_treshold"])

        valence_fractionary_correction = correction(**valence_options)
        valence_cuts, valence_gap = valence_fractionary_correction.execute()
        make_minushalf_results(valence_cuts=valence_cuts, gap=valence_gap)

    elif minushalf_yaml.correction[
            "correction_code"] == CorrectionCode.vc.name:
        valence_options["correction_indexes"] = get_simple_corrections_indexes(
            vbm_projection)
        conduction_options[
            "correction_indexes"] = get_simple_corrections_indexes(
                cbm_projection)
        valence_correction = correction(**valence_options)
        valence_cuts, _ = valence_correction.execute()
        conduction_correction = correction(**conduction_options)
        conduction_cuts, conduction_gap = conduction_correction.execute()
        make_minushalf_results(valence_cuts=valence_cuts,
                               gap=conduction_gap,
                               conduction_cuts=conduction_cuts)

    elif minushalf_yaml.correction[
            "correction_code"] == CorrectionCode.vfc.name:

        valence_options[
            "correction_indexes"] = get_fractionary_corrections_indexes(
                vbm_projection,
                treshold=minushalf_yaml.
                correction["fractionary_valence_treshold"])

        conduction_options[
            "correction_indexes"] = get_simple_corrections_indexes(
                cbm_projection)

        valence_fractionary_correction = correction(**valence_options)
        valence_cuts, _ = valence_fractionary_correction.execute()
        conduction_correction = correction(**conduction_options)
        conduction_cuts, conduction_gap = conduction_correction.execute()
        make_minushalf_results(valence_cuts=valence_cuts,
                               gap=conduction_gap,
                               conduction_cuts=conduction_cuts)

    elif minushalf_yaml.correction[
            "correction_code"] == CorrectionCode.vfcf.name:

        valence_options[
            "correction_indexes"] = get_fractionary_corrections_indexes(
                vbm_projection,
                treshold=minushalf_yaml.
                correction["fractionary_valence_treshold"])

        conduction_options[
            "correction_indexes"] = get_fractionary_corrections_indexes(
                cbm_projection,
                treshold=minushalf_yaml.
                correction["fractionary_conduction_treshold"])

        valence_fractionary_correction = correction(**valence_options)
        valence_cuts, _ = valence_fractionary_correction.execute()
        conduction_fractionary_correction = correction(**conduction_options)
        conduction_cuts, conduction_gap = conduction_fractionary_correction.execute(
        )
        make_minushalf_results(valence_cuts=valence_cuts,
                               gap=conduction_gap,
                               conduction_cuts=conduction_cuts)

    elif minushalf_yaml.correction[
            "correction_code"] == CorrectionCode.vcf.name:
        valence_options["correction_indexes"] = get_simple_corrections_indexes(
            vbm_projection)

        conduction_options[
            "correction_indexes"] = get_fractionary_corrections_indexes(
                cbm_projection,
                treshold=minushalf_yaml.
                correction["fractionary_conduction_treshold"])

        valence_correction = correction(**valence_options)
        valence_cuts, _ = valence_correction.execute()
        conduction_fractionary_correction = correction(**conduction_options)
        conduction_cuts, conduction_gap = conduction_fractionary_correction.execute(
        )
        make_minushalf_results(valence_cuts=valence_cuts,
                               gap=conduction_gap,
                               conduction_cuts=conduction_cuts)
    end_message()
