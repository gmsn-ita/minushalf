"""
Returns band-gap with the sinal changed, so
one can use minimization algorithms to find the
cut value that results in the maximum band_gap
"""
import os
import shutil
import math
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
                cut (float): Cut radius to the algorithm
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
):
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


def find_reverse_band_gap(cuts: list, *args: tuple) -> float:
    """
            Run vasp and find the eigenvalue for a value of cut

                Args:
                    cuts (float): List of cuts
                    *args (tuple): tuple containning a dictionary with the
                    following fields:
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
                    reverse_gap (float): band gap multiplied for -1

        """
    extra_args = args[0]
    cut = cuts[0]
    if cut <= extra_args["atom_potential"].vtotal.radius[1]:
        return math.inf
    runner = extra_args["runner"]
    software_factory = extra_args["software_factory"]

    find_cut_path = os.path.join(extra_args["base_path"], "find_cut")
    cut_folder = os.path.join(find_cut_path, "cut_{:.2f}".format(cut))

    if os.path.exists(cut_folder):
        shutil.rmtree(cut_folder)
    os.mkdir(cut_folder)

    software_files = extra_args["software_files"]
    for file in software_files:
        shutil.copyfile(file, os.path.join(cut_folder, file))

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

    eigenvalues = software_factory.get_eigenvalues(base_path=cut_folder)
    fermi_energy = software_factory.get_fermi_energy(base_path=cut_folder)
    atoms_map = software_factory.get_atoms_map(base_path=cut_folder)
    num_bands = software_factory.get_number_of_bands(base_path=cut_folder)
    band_projection_file = software_factory.get_band_projection_class(
        base_path=cut_folder)

    band_structure = BandStructure(eigenvalues, fermi_energy, atoms_map,
                                   num_bands, band_projection_file)

    gap_report = band_structure.band_gap()
    return (-1) * gap_report["gap"]
