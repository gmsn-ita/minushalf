"""
Returns band-gap with the sinal changed, so
one can use minimization algorithms to find the
cut value that results in the maximum band_gap
"""
import os
import re
import shutil
from loguru import logger
from subprocess import Popen, PIPE
from minushalf.io.input_file import InputFile
from minushalf.softwares.software_abstract_factory import SoftwaresAbstractFactory
from minushalf.utils.band_structure import BandStructure
from minushalf.utils.software_output import get_output_filenames


def _set_up_cut_folder(base_path: str, input_files: list, cut: float, potentials_folder: str) -> str:
    """
        Creates and populates the folder where the first principles calculations will be done

        Args:
            base_path (str): Path where the folder will be created
            cut (float): Distance to trimm the potential
            input_files (List[str]): List of input files

        Returns:
            cut_folder (str): Path to folder where the first principles calculations will be done
    """
    cut_folder = _create_cut_folder(base_path, cut)
    _copy_input_files(input_files, cut_folder, potentials_folder)
    return cut_folder


def _create_cut_folder(base_path: str, cut: float) -> str:
    """
        Creates the folder where the first principles calculations will be done

        Args:

            base_path (str): Path where the folder will be created
            cut (float): Distance to trimm the potential

        Returns:
            cut_folder (str): Path to folder where the first principles calculations will be done
    """
    find_cut_path = os.path.join(base_path, "find_cut")
    cut_folder = os.path.join(find_cut_path, "cut_{:.2f}".format(cut))

    if os.path.exists(cut_folder):
        shutil.rmtree(cut_folder)
    os.mkdir(cut_folder)

    return cut_folder


def _copy_input_files(input_files: list, destination_folder: str, potentials_folder: str) -> None:
    """
        Copy the input files for the first principles calculations

        Args:
            input_files (List[str]): List of input files
            destination_folder (str): Path to which files will be copied
    """

    for file in input_files:
        shutil.copyfile(file, os.path.join(destination_folder, file))

    # Copy previously corrected potentials
    corrected_potentials_folder = potentials_folder

    if os.path.exists(corrected_potentials_folder):
        for potential_file in os.listdir(corrected_potentials_folder):
            source = os.path.join(corrected_potentials_folder, potential_file)
            shutil.copyfile(source, os.path.join(destination_folder, potential_file))


def _get_gap(software_factory: SoftwaresAbstractFactory,
             cut_folder: str, is_indirect: bool, software_files: list) -> float:
    """
        Returns the gap value

        Args:
            software_factory (SoftwaresAbstractFactory) : Get informations for output files of the first principle calculations
            cut_folder (str): Folder where first principles calculations were made
            is_indirect (bool): This parameter determines whether a band gap calculation should be performed across various k-points.

        Returns:
            gap (float): Gap of the semiconductor material
    """

    
    filenames = get_output_filenames('QE', software_files[0])

    eigenvalues          = software_factory.get_eigenvalues(
                               filename=filenames["eigenvalues"],
                               base_path=cut_folder)
    fermi_energy         = software_factory.get_fermi_energy(
                               filename=filenames["fermi_energy"],
                               base_path=cut_folder)
    atoms_map            = software_factory.get_atoms_map(
                               filename=filenames["atoms_map"],
                               base_path=cut_folder)
    num_bands            = software_factory.get_number_of_bands(
                               filename=filenames["number_of_bands"],
                               base_path=cut_folder)
    band_projection_file = software_factory.get_band_projection_class(
                               filename=filenames["band_projection"],
                               base_path=cut_folder)

    band_structure = BandStructure(eigenvalues, fermi_energy, atoms_map,
                                   num_bands, band_projection_file)

    gap_report = band_structure.band_gap(is_indirect=is_indirect)
    return gap_report["gap"]


def _run_ld1(cwd: str, ld1_command: str) -> None:
    """
    Run ld1.x reading from INP in the given directory.
    Produces {symbol}-05.upf.tmp as output.

    Args:
        cwd (str): working directory containing INP.ldx
    """

    with open(os.path.join(cwd, "INP")) as inp:
        process = Popen(
            ld1_command,
            stdin=inp, stdout=PIPE, stderr=PIPE,
            cwd=cwd
        )
    _, stderr = process.communicate()
    if stderr:
        print(f"ld1.x stderr:\n{stderr.decode()}")

def _generate_potential(base_path: str,
                        software_factory: SoftwaresAbstractFactory,
                        potential_filename: str,
                        symbol: str,
                        orbital: str,
                        cut: float,
                        exchange_correlation_type: str,
                        max_iterations: int,
                        calculation_code: str,
                        ld1_command: str,
                        virtual_v2_command: str,
                        amplitude: float,
                        is_conduction: bool
                        ) -> None:
    """
    Routine to correct the potential using ld1.x

    Steps:
        1. Write ld1.x INP input file
        2. Run ld1.x → {symbol}-05.upf.temp
        3. Copy PP_PSWFC block from original UPF into ld1.x output
        4. Rename to the origianl potential filename
    """
    original_upf = potential_filename

    input_file = InputFile.minimum_setup(
        software="QE",
        cut=cut,
        chemical_symbol=symbol,
        exchange_correlation_code=exchange_correlation_type,
        maximum_iterations=max_iterations,
        calculation_code=calculation_code,
        file_pseudo=original_upf,
        orbital=orbital.lower(),
        amplitude=amplitude,
        is_conduction=is_conduction
    )

    inp_path = os.path.join(base_path, "INP")
    input_file.to_file(inp_path)

    _run_ld1(base_path, ld1_command)


    src = os.path.join(base_path, f"{symbol}-05.upf.temp")
    dst = os.path.join(base_path, os.path.basename(potential_filename))

    # check ld1.x output exists
    if not os.path.exists(src):
        logger.error(f"Expected ld1.x output not found: {src}")
        raise FileNotFoundError(
            f"ld1.x did not produce expected output: {src}"
        )
    
    _copy_pp_pswfc(
        source_upf=original_upf,
        target_upf=src,
    )

    shutil.move(src, dst)


def _copy_pp_pswfc(source_upf: str, target_upf: str) -> None:
    """
    Reads the PP_PSWFC block from source_upf and replaces the
    PP_PSWFC block in target_upf with it.

    pseudopotential files generated by ld1.x come with this block filled with zeros,
    sometimes impeding projwfc.x to run. 

    Args:
        source_upf (str): Path to the original unmodified UPF file.
        target_upf (str): Path to the ld1.x output UPF to be modified.
    """
    import re

    pp_pswfc_pattern = re.compile(
        r'(<PP_PSWFC[^>]*>)(.*?)(</PP_PSWFC>)',
        re.DOTALL | re.IGNORECASE
    )

    # --- read PP_PSWFC from original ---
    with open(source_upf, "r") as f:
        source_content = f.read()

    source_match = pp_pswfc_pattern.search(source_content)
    if source_match is None:
        raise ValueError(
            f"Could not find PP_PSWFC block in source UPF: {source_upf}"
        )

    # Capture the entire block including opening and closing tags
    full_pswfc_block = source_match.group(0)

    # --- read target and replace or insert ---
    with open(target_upf, "r") as f:
        target_content = f.read()

    target_match = pp_pswfc_pattern.search(target_content)

    if target_match is not None:
        # Replace existing PP_PSWFC block in target
        new_target_content = pp_pswfc_pattern.sub(
            full_pswfc_block,
            target_content,
            count=1
        )
    else:
        # No PP_PSWFC in target — insert before </UPF>
        logger.warning(
            f"No PP_PSWFC block found in target UPF: {target_upf}. "
            f"Inserting before </UPF>."
        )
        closing_upf = re.search(r'</UPF>', target_content, re.IGNORECASE)
        if closing_upf:
            insert_pos = closing_upf.start()
            new_target_content = (
                target_content[:insert_pos]
                + full_pswfc_block + "\n"
                + target_content[insert_pos:]
            )
        else:
            # No closing tag either — just append
            logger.warning(
                "No </UPF> closing tag found either. Appending PP_PSWFC at end."
            )
            new_target_content = target_content + "\n" + full_pswfc_block + "\n"

    with open(target_upf, "w") as f:
        f.write(new_target_content)

def find_negative_band_gap_qe(cuts: list, *args: tuple) -> float:
    """
                Run Qauntum ESPRESSO and return the gap value multiplied by -1

                Args:
                    cuts (float): List of cuts

                    args (tuple): tuple containning a dictionary with the fields
                                   base_path (str): Path to mkpotcar{symbol}_{orbital}
                                   symbol (str): Atom symbol
                                   default_potential_filename (str): The default potential filename for each software
                                   potfiles_folder (str): Folder containing unmodified potfiles
                                   amplitude (float): scale factor to trimming function
                                   runner (Runner): runner for the software
                                   software_factory(SoftwaresAbstractFactory): Factory for each software
                                   software_files (list): Aditional files besides potential file to make ab initio calculations

                Returns:

                    negative_gap (float): band gap multiplied for -1
    """
    extra_args = args[0]
    cut = cuts[0]
    runner = extra_args["runner"]
    ld1_command = extra_args["ld1_command"]
    virtual_v2_command = extra_args["virtual_v2_command"]
    software_factory = extra_args["software_factory"]
    is_conduction = extra_args["is_conduction"] 
    potentials_folder = os.path.join(os.path.dirname(extra_args["hidden_folder"]), "corrected_potentials")

    cut_folder = _set_up_cut_folder(extra_args["base_path"],
                                extra_args["software_files"], cut, potentials_folder) 

    _generate_potential(base_path=cut_folder,
                        potential_filename=extra_args["default_potential_filename"],
                        symbol=extra_args["symbol"],
                        orbital=extra_args["orbital"],
                        cut=cut,
                        exchange_correlation_type=extra_args["exchange_correlation_type"],
                        max_iterations=extra_args["max_iterations"],
                        calculation_code=extra_args["calculation_code"],
                        ld1_command = ld1_command,
                        virtual_v2_command=virtual_v2_command,
                        software_factory=software_factory,
                        amplitude = extra_args["amplitude"],
                        is_conduction = is_conduction)

    runner.run(cut_folder)

    is_indirect = extra_args["indirect"]
    gap = _get_gap(software_factory, cut_folder, is_indirect, extra_args["software_files"])

    # Logger
    if is_conduction:
        logger.info("CONDUCTION CORRECTION: Element {} - Orbital {}".format(
            extra_args["symbol"], extra_args["orbital"]))
        logger.info(
            "CONDUCTION CORRECTION: Current CUT value is {:.2f} a.u".format(
                cut))
        logger.info(
            "CONDUCTION CORRECTION: Current Gap value is {:.2f} eV".format(
                gap))
    else:
        logger.info("VALENCE CORRECTION: Element {} - Orbital {}".format(
            extra_args["symbol"], extra_args["orbital"]))
        logger.info(
            "VALENCE CORRECTION: Current CUT value is {:.2f} a.u".format(cut))
        logger.info(
            "VALENCE CORRECTION: Current Gap value is {:.2f} eV".format(gap))

    if cut < 0.5:
        logger.warning(
            "Pay attention, CUT values less than 0.5 a.u might have no physical meaning."
        )

    return (-1) * gap


