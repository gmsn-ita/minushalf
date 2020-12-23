"""
Implements the algorithm that automates the process
of simple valence correction and optimizes the necessary parameters.
"""
import os
import shutil
from subprocess import Popen, PIPE
import pandas as pd
from minushalf.atomic import InputFile, Vtotal
from minushalf.utils import ternary_search


class VaspValenceCorrection():
    """
    An algorithm that performs the potential correction only
    on the atom that most contributes to the character of the last valence band.
    """
    def __init__(
        self,
        root_folder: str,
        atoms: list,
        potential_filename: str,
        potential_folder: str,
        vbm_projection: pd.DataFrame,
        software_factory: any,
        runner: any,
        exchange_correlation_type: str,
        max_iterations: int,
        calculation_code: str,
        amplitude: float,
    ):
        """
        init method for the valence correction class
            Args:
                root_folder (str): Path to the folder where the valence correction will be made for each atom

                atoms (list): Atoms name

                potential_filename (str): Name of the potential file used by each
                                          software that performs ab initio calculations

                vbm_projection (pd.DataFrame): Shows the contribution of each atom to the maximum valence band

                runner: class to execute the program that makes ab initio calculations
        """
        self.root_folder = root_folder
        self.atoms = atoms
        self.potential_filename = potential_filename
        self.vbm_projection = vbm_projection
        self.potential_folder = potential_folder
        self.exchange_correlation_type = exchange_correlation_type
        self.max_iterations = max_iterations
        self.calculation_code = calculation_code
        self.runner = runner
        self.software_factory = software_factory
        self.amplitude = amplitude
        self.atom_index = None
        self.atom_potential = None

    def execute(self) -> list:
        """
        Execute valence correction algorithm
        """
        self.atom_index = self._get_correction_atom_index()
        self._generate_atom_pseudopotential()
        self._generate_occupation_potential()
        self.atom_potential = self._get_atom_potential()
        cut, eigenvalue = self._find_cut()
        return [(cut, eigenvalue)]

    def _get_atom_potential(self):
        """
        Creates atom_potential class
        """
        atom_symbol = self.atom_index[0]
        atom_folder = os.path.join(self.root_folder, atom_symbol)
        vtotal_path = os.path.join(atom_folder, "VTOTAL.ae")
        vtotal_occ_path = os.path.join(atom_folder, "VTOTAL_OCC")
        vtotal = Vtotal.from_file(vtotal_path)
        vtotal_occ = Vtotal.from_file(vtotal_occ_path)
        input_file = InputFile.minimum_setup(
            atom_symbol,
            self.exchange_correlation_type,
            self.max_iterations,
            self.calculation_code,
        )
        potential_filename = "{}.{}".format(self.potential_filename.upper,
                                            atom_symbol.lower())
        return self.software_factory.atomic_potential(
            vtotal,
            vtotal_occ,
            input_file,
            potcar_path=potential_filename,
            base_path=self.potential_folder,
        )

    def _find_cut(self) -> tuple:
        """
        Find the cut which gives the maximum gap
        """
        folder = os.path.join(self.root_folder, "find_cut")
        if os.path.exists(folder):
            shutil.rmtree(folder)
        os.mkdir(folder)
        cut, eigenvalue = ternary_search(0, 15, self._find_eigenvalue)
        return cut, eigenvalue

    def _find_eigenvalue(self, cut: float) -> float:
        """
        Run vasp and find the eigenvalue for a value of cut
            Args:
                cut (float): cut radius for the atom
            Returns:
                eigenvalue (float): eigenvalue for that specific cut
        """
        find_cut_path = os.path.join(self.root_folder, "find_cut")
        cut_folder = os.path.join(find_cut_path, "cut_{:.2f}".format(cut))

        if os.path.exists(cut_folder):
            shutil.rmtree(cut_folder)
        os.mkdir(cut_folder)

        vasp_files = ["INCAR", "KPOINTS", "POSCAR"]
        for file in vasp_files:
            shutil.copyfile(file, os.path.join(cut_folder, file))
        try:
            potential_file = open("POTCAR", "w")
            for atom in self.atoms:
                potential_filename = "{}.{}".format(
                    self.potential_filename.upper, atom.lower())
                if atom.lower() != self.atom_index[0].lower():
                    with open(potential_filename, "r") as file:
                        potential_file.write(file.read())
                else:
                    potential = self.atom_potential.correct_potential(
                        cut, self.amplitude)
                    lines = self.atom_potential.get_corrected_file_lines(
                        potential)
                    potential_file.writelines(lines)
        finally:
            potential_file.close()

        self.runner.run(cut_folder)
        band_structure = self.software_factory.band_structure(
            base_path=cut_folder)
        gap_report = band_structure.band_gap()
        return gap_report["gap"]

    def _generate_occupation_potential(self) -> None:
        """
        Generate the pseudo potential for the occupation
        of minus one half electron in the most relevant orbital
        on the VBM projection
        """
        atom_symbol = self.atom_index[0]
        secondary_quantum_number = self.atom_index[1]
        folder_path = os.path.join(self.root_folder, atom_symbol)
        if not os.path.exists(folder_path):
            raise FileNotFoundError(
                "Folder for {} does not exist".format(atom_symbol))

        process = Popen(
            ["minushalf", "occupation", "{}".format(secondary_quantum_number)],
            stdout=PIPE,
            stderr=PIPE,
            cwd=folder_path)

        _, stderr = process.communicate()
        if stderr:
            raise Exception("Call to occupation command failed")

    def _generate_atom_pseudopotential(self) -> None:
        """
        Make a dir with the atoms name,generate
        the input file and run the atomic program
        """
        atom_symbol = self.atom_index[0]
        folder_path = os.path.join(self.root_folder, atom_symbol)
        if os.path.exists(folder_path):
            shutil.rmtree(folder_path)
        input_file = InputFile.minimum_setup(
            atom_symbol,
            self.exchange_correlation_type,
            self.max_iterations,
            self.calculation_code,
        )
        input_file.to_file(folder_path)
        process = Popen(['minushalf', 'run-atomic'],
                        stdout=PIPE,
                        stderr=PIPE,
                        cwd=folder_path)

        _, stderr = process.communicate()
        if stderr:
            raise Exception('Call to atomic program failed')

    def _get_correction_atom_index(self) -> tuple:
        """
        Get dataframe index of the orbital which contributes more
        to VBM
            Returns:
                atom_index (tuple): symbol of the atom and orbital that
                contributes the most to the last band of the valence
        """
        max_projection = 0
        max_line_index = None
        max_column_index = None
        for (column_name, column_data) in self.vbm_projection.iteritems():
            current_max_projection = column_data.max()

            if current_max_projection > max_projection:
                max_projection = current_max_projection
                max_column_index = column_name
                max_line_index = column_data.idxmax()

        if not max_line_index or not max_column_index:
            raise ValueError('Projection VBM is incomplete')
        return (max_line_index, max_column_index)
