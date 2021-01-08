"""
Implements the algorithm that automates the process
simple valence correction and optimizes the necessary parameters.
"""
import itertools
import os
import shutil
from subprocess import Popen, PIPE
from scipy.optimize import minimize
import numpy as np
import pandas as pd
from minushalf.data import (CorrectionDefaultParams,
                            AtomicProgramDefaultParams, OrbitalType)
from minushalf.utils import (InputFile, Vtotal, MinushalfYaml, AtomicPotential,
                             BandStructure)
from minushalf.interfaces import (Correction, Runner, SoftwaresAbstractFactory)


class VaspValenceCorrection(Correction):
    """
    An algorithm that correct the potential in the most
    influent orbital in vbm character
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
            CorrectionDefaultParams.amplitude.name]

        self.runner = runner

        self.software_factory = software_factory

        self.corrections_index = None

        self.atom_potential = None

        self.sum_correction_percentual = 100

        self.valence_potfiles_folder = "valence_potfiles"

    def execute(self) -> tuple:
        """
        Execute valence correction algorithm
        """
        ## make valence folder in mkpotfiles
        self.root_folder = os.path.join(self.root_folder, "valence")
        if os.path.exists(self.root_folder):
            shutil.rmtree(self.root_folder)
        os.mkdir(self.root_folder)

        cuts_per_atom_orbital = {}
        self._make_valence_potential_folder()
        self.corrections_index = self._get_corrections_indexes()
        self.sum_correction_percentual = self._get_sum_correction_percentual()

        for symbol, orbitals in self.corrections_index.items():
            for orbital in orbitals:
                cut = self._find_best_correction(symbol, orbital)
                cuts_per_atom_orbital[(symbol, orbital)] = cut

        return (cuts_per_atom_orbital)

    def _make_valence_potential_folder(self):
        """
        create the folder that will store
        the potfiles corrected with an cut
        value that returns tha maximum gap
        """
        name = self.valence_potfiles_folder
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
                new_potential_file.write(potential_file.read())

            finally:
                potential_file.close()
                new_potential_file.close()

    def _get_corrections_indexes(self) -> dict:
        """
        Get index of the orbital which contributes more
        to VBM.

            Returns:
                correction_indexes (dict): A dict wherw the keys are the atoms
                symbols and the value is a list with the orbitals type to be
                corrected.
                Ex:
                {
                    'Ga': ['p','s'],
                    'N' : ['d','f'],
                }
        """

        columns_name = self.vbm_projection.columns.tolist()
        rows_name = self.vbm_projection.index.tolist()

        maximum_projection = -9999999999999999999999

        for row, column in itertools.product(rows_name, columns_name):
            projection = self.vbm_projection[column][row]
            maximum_projection = max(projection, maximum_projection)
            if np.isclose(projection, maximum_projection):
                max_row = row
                max_column = column

        correction_indexes = {max_row: [max_column]}
        return correction_indexes

    def _get_sum_correction_percentual(self) -> float:
        """
        Sum of the percentuals of orbitals
        that will be corrected.
        """
        total_sum = 0
        for symbol, orbitals in self.corrections_index.items():
            for orbital in orbitals:
                total_sum += self.vbm_projection[orbital][symbol]
        return total_sum

    def _find_best_correction(self, symbol: str, orbital: str) -> float:
        """
        Correct the potcar of the atom symbol in the
        orbital given. Then, find the best cut to the
        the correciton.

            Args:
                symbol (str): Atom symbol
                orbital (str): Orbital type (s,p,d,f)
            Returns:
                gap_and_cut (tuple): Tuple containing
                the optimum cut and the gap generated
                by the correction.
        """
        folder_name = "mkpotcar_{}_{}".format(symbol.lower(), orbital.lower())
        path = os.path.join(self.root_folder, folder_name)
        if os.path.exists(path):
            shutil.rmtree(path)
        os.mkdir(path)
        self._generate_atom_pseudopotential(path, symbol)

        percentual = 100 * (self.vbm_projection[orbital][symbol] /
                            self.sum_correction_percentual)
        self._generate_occupation_potential(path, orbital, round(percentual))
        self.atom_potential = self._get_atom_potential(path, symbol)

        cut = self._find_cut(symbol, path)
        self._write_result_in_potfile(symbol, cut, self.amplitude)
        return cut

    def _write_result_in_potfile(
        self,
        symbol: str,
        cut: float,
        amplitude: float,
    ):
        """
        Write the result of the optimization in the
        valence potfiles folder.

            Args:
                symbol (str): Atom symbol
                amplitude (float): Scale factor to trimming function
                cut (float): Cut radius to the algorithm

        """
        potential = self.atom_potential.correct_potential(cut, amplitude)
        lines = self.atom_potential.get_corrected_file_lines(potential)

        new_potential_path = os.path.join(
            self.valence_potfiles_folder,
            "{}.{}".format(self.potential_filename.upper(), symbol.lower()))

        if os.path.exists(new_potential_path):
            os.remove(new_potential_path)

        with open(new_potential_path, "w") as file:
            file.writelines(lines)

    def _generate_atom_pseudopotential(
        self,
        base_path: str,
        symbol: str,
    ) -> None:
        """
        Make a dir with the atoms name,generate
        the input file and run the atomic program

            Args:
                symbol (str): Atom symbol
                base_path (str): Path to mkpotcar{symbol}_{orbital}
        """

        folder_path = os.path.join(base_path, "pseudopotential")
        if os.path.exists(folder_path):
            shutil.rmtree(folder_path)
        os.mkdir(folder_path)
        input_file = InputFile.minimum_setup(
            symbol,
            self.exchange_correlation_type,
            self.max_iterations,
            self.calculation_code,
        )
        input_file.to_file(os.path.join(folder_path, "INP"))
        process = Popen(['minushalf', 'run-atomic', "--quiet"],
                        stdout=PIPE,
                        stderr=PIPE,
                        cwd=folder_path)

        _, stderr = process.communicate()
        if stderr:

            raise Exception("Call to atomic program failed")

    def _generate_occupation_potential(
        self,
        base_path: str,
        orbital: str,
        percentual: int,
    ) -> None:
        """
        Generate the pseudo potential for the occupation
        of a fraction of half electron.

            Args:
                base_path (str): Path to mkpotcar{symbol}_{orbital}
                orbital (str): Orbital that will be corrected
                percentual (int): Percentual of the correction
        """
        folder_path = os.path.join(base_path, "pseudopotential")
        if not os.path.exists(folder_path):

            raise FileNotFoundError(
                "Folder for pseudopotential does not exist")
        secondary_quantum_number = OrbitalType[orbital].value
        process = Popen([
            "minushalf", "occupation", "{}".format(secondary_quantum_number),
            "{}".format(percentual), "--quiet"
        ],
                        stdout=PIPE,
                        stderr=PIPE,
                        cwd=folder_path)

        _, stderr = process.communicate()
        if stderr:
            raise Exception("Call to occupation command failed")

    def _get_atom_potential(self, base_path: str,
                            symbol: str) -> AtomicPotential:
        """
        Creates atom_potential class

            Args:
                symbol (str): Atom symbol
                base_path (str): Path to mkpotcar{symbol}_{orbital}
        """
        pseudopotential_folder = os.path.join(base_path, "pseudopotential")
        vtotal_path = os.path.join(pseudopotential_folder, "VTOTAL.ae")
        vtotal_occ_path = os.path.join(pseudopotential_folder, "VTOTAL_OCC")
        vtotal = Vtotal.from_file(vtotal_path)
        vtotal_occ = Vtotal.from_file(vtotal_occ_path)

        potential_filename = "{}.{}".format(self.potential_filename.upper(),
                                            symbol.lower())

        potential_class = self.software_factory.get_potential_class(
            potential_filename, self.valence_potfiles_folder)

        return AtomicPotential(vtotal, vtotal_occ, potential_class)

    def _find_cut(self, symbol: str, base_path: str) -> tuple:
        """
        Find the cut which gives the maximum gap

            Args:
                symbol (str): Atom symbol
                base_path (str): Path to mkpotcar{symbol}_{orbital}
        """
        folder = os.path.join(base_path, "find_cut")
        if os.path.exists(folder):
            shutil.rmtree(folder)
        os.mkdir(folder)
        function_args = {
            "base_path": base_path,
            "software_factory": self.software_factory,
            "runner": self.runner,
            "symbol": symbol,
            "default_potential_filename": self.potential_filename,
            "atom_potential": self.atom_potential,
            "potfiles_folder": self.potential_folder,
            "amplitude": self.amplitude,
            "atoms": self.atoms,
        }
        res = minimize(VaspValenceCorrection.find_band_gap,
                       x0=3,
                       args=(function_args),
                       method="Nelder-Mead",
                       tol=0.01)
        cut = res.x[0]

        if not res.success:
            print("Optimization failed")
            raise Exception("Optimization failed.")

        return cut

    @staticmethod
    def _get_corrected_potfile_lines(
        cut: float,
        atom_potential: AtomicPotential,
        amplitude: float,
    ) -> list:
        """
        Modify the potential file that was corrected

            Args:
                cut (float): Cut radius to the algorithm
                atom_potential(AtomicPotential): Class that represents the
                potential of the atom
                amplitude (float): Scale factor to trimming function
        """
        potential = atom_potential.correct_potential(cut, amplitude)
        lines = atom_potential.get_corrected_file_lines(potential)
        return lines

    @staticmethod
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

    @staticmethod
    def find_band_gap(cut: float, *args: tuple) -> float:
        """
            Run vasp and find the eigenvalue for a value of cut

                Args:
                    cut (float): Cut radius to the algorithm
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

                Returns:
                    reverse_gap (float): band gap multiplied for -1

        """
        extra_args = args[0]
        runner = extra_args["runner"]
        software_factory = extra_args["software_factory"]

        find_cut_path = os.path.join(extra_args["base_path"], "find_cut")
        cut_folder = os.path.join(find_cut_path, "cut_{:.2f}".format(cut))

        if os.path.exists(cut_folder):
            shutil.rmtree(cut_folder)
        os.mkdir(cut_folder)

        vasp_files = ["INCAR", "KPOINTS", "POSCAR"]
        for file in vasp_files:
            shutil.copyfile(file, os.path.join(cut_folder, file))

        potfile_lines = VaspValenceCorrection._get_corrected_potfile_lines(
            cut,
            extra_args["atom_potential"],
            extra_args["amplitude"],
        )
        VaspValenceCorrection._join_potfiles(
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
