"""
Implements the algorithm that automates the process
of vasp correction and optimizes the necessary parameters.
"""
import os
import shutil
from subprocess import Popen, PIPE
from scipy.optimize import minimize
from loguru import logger
import pandas as pd
from minushalf.data import (CorrectionDefaultParams,
                            AtomicProgramDefaultParams, OrbitalType)
from minushalf.utils import (InputFile, Vtotal, MinushalfYaml, AtomicPotential,
                             BandStructure, find_reverse_band_gap)
from minushalf.interfaces import (Correction, Runner, SoftwaresAbstractFactory)


class VaspCorrection(Correction):
    """
    An algorithm that realizes corrections for
    VASP software
    """
    def __init__(
        self,
        root_folder: str,
        software_factory: SoftwaresAbstractFactory,
        runner: Runner,
        minushalf_yaml: MinushalfYaml,
        band_projection: pd.DataFrame,
        atoms: list,
        is_conduction: bool,
        correction_indexes: dict,
    ):
        """
        init method for the vasp correction class
            Args:
                root_folder (str): Path to the folder where the  correction will be made for each atom

                atoms (list): Atoms name

                potential_filename (str): Name of the potential file used by each
                                          software that performs ab initio calculations

                band_projection (pd.DataFrame): Shows the contribution of each atom in the CBM or VBM

                runner: class to execute the program that makes ab initio calculations
        """
        self.root_folder = root_folder

        self.atoms = atoms

        self.potential_filename = software_factory.get_potential_class(
        ).get_name()

        self.band_projection = band_projection

        if is_conduction:
            self.potential_folder = "corrected_valence_potfiles"
        else:
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

        if is_conduction:
            self.cut_initial_guess = minushalf_yaml.correction[
                CorrectionDefaultParams.conduction_cut_guess.name]
        else:
            self.cut_initial_guess = minushalf_yaml.correction[
                CorrectionDefaultParams.valence_cut_guess.name]

        self.tolerance = minushalf_yaml.correction[
            CorrectionDefaultParams.tolerance.name]

        self.runner = runner

        self.software_factory = software_factory

        self.atom_potential = None

        self.sum_correction_percentual = 100

        if is_conduction:
            self.corrected_potfiles_folder = "corrected_conduction_potfiles"
        else:
            self.corrected_potfiles_folder = "corrected_valence_potfiles"

        if is_conduction:
            self.correction_type = "conduction"
        else:
            self.correction_type = "valence"

        self.is_conduction = is_conduction

        self.correction_indexes = correction_indexes

        self.software_files = ["INCAR", "POSCAR", "KPOINTS"]

    @property
    def potential_folder(self) -> str:
        """
        Returns:
            Name of the folder that helds all the potential files
            initially not corrected.
        """
        return self._potential_folder

    @potential_folder.setter
    def potential_folder(self, path: str) -> None:
        """
        Verify if the folder exists and contains all files needed

        Args:
            path (str): Path of the folder that helds all the potential files
            initially not corrected.
        """

        for atom in self.atoms:
            filename = "{}.{}".format(self.potential_filename.upper(),
                                      atom.lower())
            abs_path = os.path.join(path, filename)
            if not os.path.exists(abs_path):
                logger.error("Potential folder incomplete")
                raise FileNotFoundError("Potential folder lacks of files.")

        self._potential_folder = path

    def execute(self) -> tuple:
        """
        Execute vasp correction algorithm
        """
        ## make (valence|conduction) folder in mkpotfiles
        self.root_folder = os.path.join(self.root_folder, self.correction_type)
        if os.path.exists(self.root_folder):
            shutil.rmtree(self.root_folder)
        os.mkdir(self.root_folder)

        cuts_per_atom_orbital = {}
        self._make_corrected_potential_folder()

        self.sum_correction_percentual = self._get_sum_correction_percentual()

        for symbol, orbitals in self.correction_indexes.items():
            for orbital in orbitals:
                cut = self._find_best_correction(symbol, orbital)
                cuts_per_atom_orbital[(symbol, orbital)] = cut
        gap = self._get_result_gap()
        return (cuts_per_atom_orbital, gap)

    def _get_result_gap(self) -> float:
        """
        Return the gap after the optimization of all potfiles
        """
        calculation_folder = "calculate_{}_gap".format(self.correction_type)
        if os.path.exists(calculation_folder):
            shutil.rmtree(calculation_folder)
        os.mkdir(calculation_folder)

        for file in self.software_files:
            shutil.copyfile(file, os.path.join(calculation_folder, file))

        potfile_path = os.path.join(calculation_folder,
                                    self.potential_filename)
        potential_file = open(potfile_path, "w")
        try:
            for atom in self.atoms:
                atom_potfilename = "{}.{}".format(
                    self.potential_filename.upper(), atom.lower())
                atom_potpath = os.path.join(self.corrected_potfiles_folder,
                                            atom_potfilename)
                with open(atom_potpath) as file:
                    potential_file.write(file.read())
        finally:
            potential_file.close()
        ## Run ab initio calculations
        self.runner.run(calculation_folder)

        eigenvalues = self.software_factory.get_eigenvalues(
            base_path=calculation_folder)
        fermi_energy = self.software_factory.get_fermi_energy(
            base_path=calculation_folder)
        atoms_map = self.software_factory.get_atoms_map(
            base_path=calculation_folder)
        num_bands = self.software_factory.get_number_of_bands(
            base_path=calculation_folder)
        band_projection_file = self.software_factory.get_band_projection_class(
            base_path=calculation_folder)

        band_structure = BandStructure(eigenvalues, fermi_energy, atoms_map,
                                       num_bands, band_projection_file)

        gap_report = band_structure.band_gap()
        return gap_report["gap"]

    def _make_corrected_potential_folder(self):
        """
        create the folder that will store
        the potfiles corrected with an cut
        value that returns tha maximum gap
        """
        name = self.corrected_potfiles_folder
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

    def _get_sum_correction_percentual(self) -> float:
        """
        Sum of the percentuals of orbitals
        that will be corrected.
        """
        total_sum = 0
        for symbol, orbitals in self.correction_indexes.items():
            for orbital in orbitals:
                total_sum += self.band_projection[orbital][symbol]
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

        percentual = 100 * (self.band_projection[orbital][symbol] /
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
        corrected potfiles folder.

            Args:
                symbol (str): Atom symbol
                amplitude (float): Scale factor to trimming function
                cut (float): Cut radius to the algorithm

        """
        potential = self.atom_potential.correct_potential(
            cut, amplitude, is_conduction=self.is_conduction)
        lines = self.atom_potential.get_corrected_file_lines(potential)

        new_potential_path = os.path.join(
            self.corrected_potfiles_folder,
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
            potential_filename, self.corrected_potfiles_folder)

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
            "software_files": self.software_files,
            "is_conduction": self.is_conduction,
        }
        res = minimize(find_reverse_band_gap,
                       x0=self.cut_initial_guess,
                       args=(function_args),
                       method="Nelder-Mead",
                       tol=self.tolerance)
        cut = res.x[0]

        if not res.success:
            print("Optimization failed")
            raise Exception("Optimization failed.")

        return cut
