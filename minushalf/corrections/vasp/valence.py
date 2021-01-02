"""
Implements the algorithm that automates the process
of simple valence correction and optimizes the necessary parameters.
"""
import os
import shutil
from subprocess import Popen, PIPE
import pandas as pd
from minushalf.data import (CorrectionDefaultParams,
                            AtomicProgramDefaultParams)
from minushalf.utils import (ternary_search, InputFile, Vtotal, MinushalfYaml,
                             AtomicPotential, BandStructure)
from minushalf.interfaces import (Correction, Runner, SoftwaresAbstractFactory)


class VaspValenceCorrection(Correction):
    """
    An algorithm that performs the potential correction only
    on the atom that most contributes to the character of the last valence band.
    """
    def __init__(
        self,
        root_folder: str,
        software_factory: SoftwaresAbstractFactory,
        runner: Runner,
        minushalf_yaml: MinushalfYaml,
        vbm_projection: pd.DataFrame,
        atoms: list,
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

        self.potential_filename = software_factory.get_potential_class(
        ).get_name()

        self.vbm_projection = vbm_projection

        self.potential_folder = minushalf_yaml.correction[
            CorrectionDefaultParams.potfiles_folder.name]

        self.exchange_correlation_type = minushalf_yaml.atomic_program[
            AtomicProgramDefaultParams.exchange_correlation_code.name]

        self.max_iterations = minushalf_yaml.atomic_program[
            AtomicProgramDefaultParams.max_iterations.name]

        self.calculation_code = minushalf_yaml.atomic_program[
            AtomicProgramDefaultParams.calculation_code.name]

        self.amplitude = minushalf_yaml.correction[
            CorrectionDefaultParams.amplitude.value]

        self.runner = runner

        self.software_factory = software_factory

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

        potential_filename = "{}.{}".format(self.potential_filename.upper,
                                            atom_symbol.lower())
        potential_class = self.software_factory.get_potential_class(
            potential_filename, self.potential_folder)

        return AtomicPotential(vtotal, vtotal_occ, potential_class)

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

        self._make_potcar_folder(cut)
        self._join_valence_potfiles(cut_folder)

        self.runner.run(cut_folder)

        eigenvalues = self.software_factory.get_eigenvalues(
            base_path=cut_folder)
        fermi_energy = self.software_factory.get_fermi_energy(
            base_path=cut_folder)
        atoms_map = self.software_factory.get_atoms_map(base_path=cut_folder)
        num_bands = self.software_factory.get_number_of_bands(
            base_path=cut_folder)
        band_projection_file = self.software_factory.get_band_projection_class(
            base_path=cut_folder)

        band_structure = BandStructure(eigenvalues, fermi_energy, atoms_map,
                                       num_bands, band_projection_file)

        gap_report = band_structure.band_gap()
        return gap_report["gap"]

    def _make_potcar_folder(
        self,
        cut: float,
        name: str = "valence_potfiles",
    ):
        """
        Make potcar folder for each cut calculation
        """
        if os.path.exists(name):
            shutil.rmtree(name)
        os.mkdir(name)
        for atom in self.atoms:
            try:
                potential_filename = "{}.{}".format(
                    self.potential_filename.upper(), atom.lower())
                potential_path = os.path.join(self.potential_folder,
                                              potential_filename)
                potential_file = open(potential_path, "r")

                new_potential_path = os.path.join(name, potential_filename)
                new_potential_file = open(new_potential_path, "w")

                if atom.lower() != self.atom_index[0].lower():

                    new_potential_file.write(potential_file.read())
                else:
                    potential = self.atom_potential.correct_potential(
                        cut, self.amplitude)
                    lines = self.atom_potential.get_corrected_file_lines(
                        potential)
                    new_potential_file.writelines(lines)
            finally:
                potential_file.close()
                new_potential_file.close()

    def _join_valence_potfiles(self,
                               base_path: str,
                               pot_folder: str = "valence_potfiles"):
        """
        Join valence potfiles in one
        """
        path = os.path.join(base_path, self.potential_filename.upper())
        if os.path.exists(path):
            shutil.rmtree(path)
        os.mkdir(path)
        potcar_file = open(path, "w")

        try:
            for atom in self.atoms:

                potential_filename = "{}.{}".format(
                    self.potential_filename.upper(), atom.lower())
                potential_path = os.path.join(pot_folder, potential_filename)

                with open(potential_path, "r") as file:
                    potcar_file.write(file.read())
        finally:
            potcar_file.close()
