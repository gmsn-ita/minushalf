"""
Returns band-gap with the sinal changed, so
one can use minimization algorithms to find the
cut value that results in the maximum band_gap
"""
import os
import shutil
import math
from loguru import logger
from minushalf.interfaces import SoftwaresAbstractFactory
from .atomic_potential import AtomicPotential
from .band_structure import BandStructure


def _get_corrected_potfile_lines(
    cut: float,
    atom_potential: AtomicPotential,
    amplitude: float,
    is_conduction: bool,
) -> list:
    """
        Modify the potential file that was corrected

            Args:

                cut (float): Distance to trimm the potential
                atom_potential(AtomicPotential): Class that represents the
                                                 potential of the atom
                amplitude (float): Scale factor to trimming function
        """
    potential = atom_potential.correct_potential(cut,
                                                 amplitude,
                                                 is_conduction=is_conduction)
    lines = atom_potential.get_corrected_file_lines(potential)
    return lines


def _join_potfiles(
    base_path: str,
    default_potential_filename: str,
    atoms: list,
    potfiles_folder: str,
    symbol: str,
    corrected_potfile_lines: list,
) -> None:
    """
        Join valence potfiles in one

        Args:

            base_path (str): Path to mkpotcar{symbol}_{orbital}
            symbol (str): Atom symbol
            default_potential_filename (str): The default potential filename for each software
            potfiles_folder (str): Folder containing unmodified potfiles
            symbol (str): symbol of the atom to be corrected
            potfile_lines (list): Lines of the corrected potential file

        """
    path = os.path.join(base_path, default_potential_filename.upper())
    if os.path.exists(path):
        os.remove(path)
    potential_file = open(path, "w")

    try:
        for atom in atoms:

            if not atom == symbol:
                potential_filename = "{}.{}".format(
                    default_potential_filename.upper(), atom.lower())
                potential_path = os.path.join(potfiles_folder,
                                              potential_filename)

                with open(potential_path, "r") as file:
                    potential_file.write(file.read())
            else:
                potential_file.writelines(corrected_potfile_lines)
    finally:
        potential_file.close()


def _set_up_cut_folder(base_path: str, input_files: list, cut: float) -> str:
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
    _copy_input_files(input_files, cut_folder)
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


def _copy_input_files(input_files: list, destination_folder: str) -> None:
    """
        Copy the input files for the first principles calculations

        Args:
            input_files (List[str]): List of input files
            destination_folder (str): Path to which files will be copied
    """

    for file in input_files:
        shutil.copyfile(file, os.path.join(destination_folder, file))


def _get_gap(software_factory: SoftwaresAbstractFactory,
             cut_folder: str) -> float:
    """
        Returns the gap value

        Args:
            software_factory (SoftwaresAbstractFactory) : Get informations for output files of the first principle calculations
            cut_folder (str): Folder where first principles calculations were made
        
        Returns:
            gap (float): Gap of the semiconductor material
    """
    eigenvalues = software_factory.get_eigenvalues(base_path=cut_folder)
    fermi_energy = software_factory.get_fermi_energy(base_path=cut_folder)
    atoms_map = software_factory.get_atoms_map(base_path=cut_folder)
    num_bands = software_factory.get_number_of_bands(base_path=cut_folder)
    band_projection_file = software_factory.get_band_projection_class(
        base_path=cut_folder)

    band_structure = BandStructure(eigenvalues, fermi_energy, atoms_map,
                                   num_bands, band_projection_file)

    gap_report = band_structure.band_gap()
    return gap_report["gap"]


def find_negative_band_gap(cuts: list, *args: tuple) -> float:
    """
                Run vasp and return the gap value multiplied by -1

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
                                   atom_potential(AtomicPotential): Holds fourier transforms of the potential
                                   software_files (list): Aditional files besides potential file to make ab initio calculations

                Returns:

                    negative_gap (float): band gap multiplied for -1
    """
    extra_args = args[0]
    cut = cuts[0]
    runner = extra_args["runner"]
    software_factory = extra_args["software_factory"]
    inplace = extra_args["inplace"]
    is_conduction = extra_args["is_conduction"]

    if cut <= extra_args["atom_potential"].vtotal.radius[1]:
        return math.inf

    ## Add condition to include inplace calculations
    if not inplace:
        cut_folder = _set_up_cut_folder(extra_args["base_path"],
                                        extra_args["software_files"], cut)
    else:
        cut_folder = '.'

    potfile_lines = _get_corrected_potfile_lines(
        cut,
        extra_args["atom_potential"],
        extra_args["amplitude"],
        extra_args["is_conduction"],
    )
    _join_potfiles(
        cut_folder,
        extra_args["default_potential_filename"],
        extra_args["atoms"],
        extra_args["potfiles_folder"],
        extra_args["symbol"],
        potfile_lines,
    )

    runner.run(cut_folder)

    gap = _get_gap(software_factory, cut_folder)

    ## Logger
    if is_conduction:
        logger.info(
            "CONDUCTION CORRECTION: Current CUT value is {:.2f} a.u".format(
                cut))
        logger.info(
            "CONDUCTION CORRECTION: Current Gap value is {:.2f} eV".format(
                gap))
    else:
        logger.info(
            "VALENCE CORRECTION: Current CUT value is {:.2f} a.u".format(cut))
        logger.info(
            "VALENCE CORRECTION: Current Gap value is {:.2f} eV".format(gap))

    if cut < 0.5:
        logger.warning(
            "Pay attention, CUT values less than 0.5 a.u might have no physical meaning."
        )

    return (-1) * gap
