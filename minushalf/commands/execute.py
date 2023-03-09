"""
Execute command
"""
import os
import sys
import shutil
from collections import OrderedDict
from typing import List
import numpy as np
import click
from loguru import logger
from minushalf.io.minushalf_yaml import MinushalfYaml
from minushalf.io.make_minushalf_results import make_minushalf_results
from minushalf.utils.cli_messages import welcome_message,end_message
from minushalf.utils.get_correction_params import get_valence_correction_params, get_conduction_correction_params

from minushalf.softwares.vasp import (Vasp)
from minushalf.corrections.correction import (DFTCorrection)
from minushalf.data.minushalf_yaml_default_configuration import CorrectionDefaultParams
from minushalf.data.softwares import (Softwares)

from minushalf.interfaces.software_abstract_factory import (SoftwaresAbstractFactory)


def get_atoms_list(factory: SoftwaresAbstractFactory) -> list:
    """
    Returns atoms_list
    """
    atoms_map = factory.get_atoms_map()
    atoms = [atoms_map[key] for key in sorted(atoms_map)]
    return list(OrderedDict.fromkeys(atoms))


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
    correction_factory_chooser = {Softwares.vasp.value: DFTCorrection}
    software_factory_chooser = {Softwares.vasp.value: Vasp()}

    software_name = minushalf_yaml.get_software_name()
    correction = correction_factory_chooser[software_name]
    software_factory = software_factory_chooser[software_name]

    ## Makes abinition calculation
    logger.info("Running ab initio calculations")
    software_configurations = minushalf_yaml.get_software_configurations_params(
    )
    runner = software_factory.get_runner(**software_configurations)
    runner.run()

    ## Makes root folder
    logger.info("Make potfiles folder")
    root_folder = "mkpotfiles"
    if os.path.exists(root_folder):
        shutil.rmtree(root_folder)
    os.mkdir(root_folder)

    ## get atoms list
    logger.info("Get atoms list")
    atoms = get_atoms_list(software_factory)

    ## amplitude logger
    if not np.isclose(minushalf_yaml.get_amplitude(),
                      CorrectionDefaultParams.amplitude.value):
        logger.warning(
            "Amplitude value is different from 1.0. This is not recommended unless you know exactly what you are doing."
        )

    valence_options = get_valence_correction_params(minushalf_yaml,
                                                    software_factory,
                                                    atoms=atoms,
                                                    runner=runner,
                                                    root_folder=root_folder)

    conduction_options = get_conduction_correction_params(
        minushalf_yaml,
        software_factory,
        atoms=atoms,
        runner=runner,
        root_folder=root_folder)

    correction_code = minushalf_yaml.get_correction_code()

    logger.info("Doing corrections")
    gap = None
    valence_cuts = None
    conduction_cuts = None

    if "v" in correction_code:
        valence_correction = correction(**valence_options)
        valence_cuts, valence_gap = valence_correction.execute()
        gap = valence_gap
        make_minushalf_results(valence_cuts=valence_cuts, gap=valence_gap)

    if "c" in correction_code:
        conduction_correction = correction(**conduction_options)
        conduction_cuts, conduction_gap = conduction_correction.execute()
        gap = conduction_gap

    make_minushalf_results(valence_cuts=valence_cuts,
                           gap=gap,
                           conduction_cuts=conduction_cuts)

    end_message()
