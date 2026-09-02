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
from minushalf.utils.orbital import (OrbitalType)
from minushalf.io.vtotal import Vtotal
from minushalf.io.input_file import InputFile
from minushalf.utils.atomic_potential import AtomicPotential
from minushalf.utils.band_structure import BandStructure
from minushalf.utils.negative_band_gap import find_negative_band_gap
from minushalf.utils.negative_band_gap_qe import find_negative_band_gap_qe
from minushalf.io.correction_interface import Correction
from minushalf.softwares.runner import Runner
from minushalf.softwares.software_abstract_factory import SoftwaresAbstractFactory
from minushalf.utils.software_output import get_output_filenames


class DFTCorrection(Correction):
    """
    An algorithm that realizes corrections for
    VASP software
    """

    def __init__(
        self,
        root_folder: str,
        potential_filename: str,
        potential_folder: str,
        exchange_correlation_type: str,
        max_iterations: int,
        software_factory: SoftwaresAbstractFactory,
        runner: Runner,
        calculation_code: str,
        amplitude: float,
        cut_initial_guess: dict,
        tolerance: float,
        input_files: list,
        indirect: bool,
        corrected_potfiles_folder: str,
        correction_type: str,
        band_projection: pd.DataFrame,
        atoms: list,
        is_conduction: bool,
        correction_indexes: dict,
        divide_character: list,
        **kwargs
    ):
        """
        init method for the vasp correction class
            Args:
                root_folder (str): Path to the folder where the  correction will be made for each atom

                atoms (list): Atoms name

                potential_filename (str): Name of the potential file used by each
                                          software that performs ab initio calculations

                band_projection (pd.DataFrame): Shows the contribution of each atom in the CBM or VBM

                runner (Runner): class to execute the program that makes ab initio calculations

                only_conduction (bool): Conduction correction without previous valence correction

                indirect (bool): Realize calculations considering indirect gaps
        """
        self.root_folder = root_folder

        self.atoms = atoms

        self.potential_filename = potential_filename

        self.band_projection = band_projection

        self.potential_folder = potential_folder

        self.exchange_correlation_type = exchange_correlation_type

        self.max_iterations = max_iterations

        self.calculation_code = calculation_code

        self.amplitude = amplitude

        self.cut_initial_guess = cut_initial_guess

        self.tolerance = tolerance

        self.runner = runner

        self.software_factory = software_factory

        self.atom_potential = None

        self.sum_correction_percentual = 100

        self.corrected_potfiles_folder = corrected_potfiles_folder

        self.correction_type = correction_type

        self.is_conduction = is_conduction

        self.correction_indexes = correction_indexes

        self.input_files = input_files

        self.indirect = indirect

        self.divide_character = divide_character

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
        # make (valence|conduction) folder in mkpotfiles
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

        gap = self._get_result_gap(is_indirect=self.indirect)
        return (cuts_per_atom_orbital, gap)

    def _get_result_gap(self, is_indirect: bool) -> float:
        """
        Return the gap after the optimization of all potfiles
        """
        calculation_folder = "./.minushalf/calculate_{}_gap".format(
            self.correction_type)
        if os.path.exists(calculation_folder):
            shutil.rmtree(calculation_folder)
        os.mkdir(calculation_folder)

        for file in self.input_files:
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
        # Run ab initio calculations
        self.runner.run(calculation_folder)

        filenames = get_output_filenames('VASP')

        eigenvalues = self.software_factory.get_eigenvalues(
            base_path=calculation_folder,
            filename=filenames["eigenvalues"])
        fermi_energy = self.software_factory.get_fermi_energy(
            base_path=calculation_folder,
            filename=filenames["fermi_energy"])
        atoms_map = self.software_factory.get_atoms_map(
            base_path=calculation_folder,
            filename=filenames["atoms_map"])
        num_bands = self.software_factory.get_number_of_bands(
            base_path=calculation_folder,
            filename=filenames["number_of_bands"])
        band_projection_file = self.software_factory.get_band_projection_class(
            base_path=calculation_folder,
            filename=filenames["band_projection"])

        band_structure = BandStructure(eigenvalues, fermi_energy, atoms_map,
                                       num_bands, band_projection_file)

        gap_report = band_structure.band_gap(is_indirect)
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

        if total_sum == 0:
            logger.error(
                "No orbital selected for correction. Check you treshhold")
            raise ValueError(
                "No orbital selected for correction. Check you treshhold")

        return total_sum

    def _find_best_correction(self, symbol: str, orbital: str) -> tuple:
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
        folder_name = f"mkpotcar_{symbol.lower()}_{orbital.lower()}"
        path = os.path.join(self.root_folder, folder_name)
        if os.path.exists(path):
            shutil.rmtree(path)
        os.mkdir(path)
        self._generate_atom_potential(path, symbol)

        percentuals = {}
        # Check for bonds with equal atoms
        number_equal_neighbors = self.divide_character[(symbol.capitalize(),
                                                        orbital.lower())]

        value = (100 / (1 + number_equal_neighbors)) * (
            self.band_projection[orbital][symbol] /
            self.sum_correction_percentual)
        logger.info(f"percentual of half eletron is {round(value)}")
        percentuals[orbital] = round(value)

        self._generate_occupation_potential(path, percentuals)

        self.atom_potential = self._get_atom_potential_class(path, symbol)

        cut = self._find_cut(symbol, path, orbital)
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

    def _generate_atom_potential(
        self,
        base_path: str,
        symbol: str,
    ) -> None:
        """
        Make a dir with the atoms name,generate
        the input file and run the atomic program

            Args:
                symbol (str): Atom symbol
                base_path (str): Path to mkpotcar{symbol}
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
        percentuals: dict,
    ) -> None:
        """
        Generate the potential for the occupation
        of a fraction of half electron.

            Args:
                base_path (str): Path to mkpotcar{symbol}
                percentual (dict): Dict where the key is the orbital type
                and the value is the percentual to be corrected
        """
        folder_path = os.path.join(base_path, "pseudopotential")
        if not os.path.exists(folder_path):

            raise FileNotFoundError(
                "Folder for pseudopotential does not exist")
        secondary_quantum_numbers = ",".join(
            [str(OrbitalType[key].value) for key in percentuals.keys()])

        joined_percentual = ",".join(
            [str(value) for value in percentuals.values()])

        process = Popen([
            "minushalf", "occupation", "{}".format(secondary_quantum_numbers),
            "{}".format(joined_percentual), "--quiet"
        ],
            stdout=PIPE,
            stderr=PIPE,
            cwd=folder_path)

        _, stderr = process.communicate()
        if stderr:
            raise Exception("Call to occupation command failed")

    def _get_atom_potential_class(self, base_path: str,
                                  symbol: str) -> AtomicPotential:
        """
        Creates atom_potential class

            Args:
                symbol (str): Atom symbol
                base_path (str): Path to mkpotcar{symbol}
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

    def _find_cut(self, symbol: str, base_path: str, orbital: str) -> tuple:
        """
        Find the cut which gives the maximum gap

            Args:
                symbol (str): Atom symbol
                base_path (str): Path to mkpotcar{symbol}
        """
        if not self.indirect:
            folder = os.path.join(base_path, "find_cut")
            if os.path.exists(folder):
                shutil.rmtree(folder)
            os.mkdir(folder)
        else:
            base_path = '.'  # Local Path

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
            "software_files": self.input_files,
            "is_conduction": self.is_conduction,
            "indirect": self.indirect,
            "orbital": orbital
        }
        cut_initial_guess = self.cut_initial_guess[(symbol.capitalize(),
                                                    orbital.lower())]

        result = minimize(find_negative_band_gap,
                          x0=cut_initial_guess,
                          args=(function_args),
                          method="Nelder-Mead",
                          options={'xatol': self.tolerance})

        if not result.success:
            logger.error("Optimization failed")
            raise Exception("Optimization failed.")

        cut = result.x[0]
        return cut


class QECorrection(Correction):
    """
    An algorithm that realizes DFT-1/2 corrections for
    Quantum ESPRESSO software.
    """

    def __init__(
        self,
        root_folder: str,
        potential_filename: str,
        potential_folder: str,         # folder with uncorrected .upf files per element
        exchange_correlation_type: str,
        max_iterations: int,
        software_factory: SoftwaresAbstractFactory,
        runner: Runner,
        calculation_code: str,
        amplitude: float,
        cut_initial_guess: dict,
        tolerance: float,
        input_files: list,             # QE input files (scf.in, etc.)
        indirect: bool,
        corrected_potfiles_folder: str,
        correction_type: str,
        band_projection: pd.DataFrame,
        atoms: list,
        is_conduction: bool,
        correction_indexes: dict,
        divide_character: list,
        commands: dict,
    ):
        """
        init method for the qe correction class
            Args:
                root_folder (str): Path to the folder where the  correction will be made for each atom

                potential_filename (str): Path of the potential choosen for correction.

                potential_folder (str): Folder containing all potential files

                atoms (list): Atoms name

                band_projection (pd.DataFrame): Shows the contribution of each atom in the CBM or VBM

                runner (Runner): class to execute the program that makes ab initio calculations

                only_conduction (bool): Conduction correction without previous valence correction

                indirect (bool): Realize calculations considering indirect gaps

                input_files (list): Input files for QE, in the following order [scf.in, A1.upf, A2.upf, ...]

                commands (list): List of commands to run pw.x, projwfc.x, ld1.x and vistual_v2.x
        """
        self.root_folder              = root_folder
        self.atoms                    = atoms
        self.potential_filename       = potential_filename
        self.band_projection          = band_projection
        self.potential_folder         = potential_folder
        self.exchange_correlation_type = exchange_correlation_type
        self.max_iterations           = max_iterations
        self.calculation_code         = calculation_code
        self.amplitude                = amplitude
        self.cut_initial_guess        = cut_initial_guess
        self.tolerance                = tolerance
        self.runner                   = runner
        self.software_factory         = software_factory
        self.atom_potential           = None
        self.sum_correction_percentual = 100
        self.corrected_potfiles_folder = corrected_potfiles_folder
        self.correction_type          = correction_type
        self.is_conduction            = is_conduction
        self.correction_indexes       = correction_indexes
        self.input_files              = input_files
        self.indirect                 = indirect
        self.divide_character         = divide_character
        self.commands                 = commands

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
        Verify that the potential file exists inside the potential folder.

        The potential filename is provided by self.potential_filename.
        """
        abs_path = os.path.join(path, self.potential_filename)

        if not os.path.exists(abs_path):
            logger.error("Potential folder incomplete")
            raise FileNotFoundError(
                f"Potential folder lacks {self.potential_filename}."
            )

        self._potential_folder = path


    def _backup_original_potentials(self) -> None:
        """
        Backup the original pseudopotential files.
        """
        main_folder = os.path.dirname(os.path.abspath(self.input_files[0]))


        backup_folder = os.path.join(
            main_folder,
            "minushalf_original_potentials"
        )

        if os.path.exists(backup_folder):
            shutil.rmtree(backup_folder)

        os.makedirs(backup_folder)

        for potential in self.input_files[1:]:
            source = os.path.abspath(potential)
            destination = os.path.join(
                backup_folder,
                os.path.basename(potential)
            )

            if not os.path.exists(source):
                raise FileNotFoundError(
                    f"Original pseudopotential not found: {source}"
                )

            shutil.copy2(source, destination)

        logger.info(
            f"Original pseudopotentials backed up to {backup_folder}"
        )


    def _copy_corrected_potential(
        self,
        symbol: str,
        orbital: str,
        correction_folder: str
    ) -> None:
        """
        Copy the corrected pseudopotential generated by ld1.x to the
        main QE calculation folder.
        """
        main_folder = os.path.dirname(os.path.abspath(self.input_files[0]))

        # Find the pseudopotential corresponding to this atom.
        potential_filename = get_output_filenames(
            "QE",
            self.input_files[0],
            atom=symbol
        )["potential"]

        corrected_potential = os.path.join(
            correction_folder,
            f"{symbol}-05.upf"
        )

        destination = os.path.join(
            main_folder,
            os.path.basename(potential_filename)
        )

        if not os.path.exists(corrected_potential):
            raise FileNotFoundError(
                f"Corrected pseudopotential not found: "
                f"{corrected_potential}"
            )

        shutil.copy2(
            corrected_potential,
            destination
        )

        logger.info(
            f"Corrected {symbol} {orbital} potential copied to "
            f"{destination}"
        )



    def execute(self) -> tuple:
        """
        Execute the Quantum ESPRESSO correction algorithm.
        
        """
        self._backup_original_potentials()

        self.root_folder = os.path.join(self.root_folder, self.correction_type)
        if os.path.exists(self.root_folder):
            shutil.rmtree(self.root_folder)
        os.mkdir(self.root_folder)

        cuts_per_atom_orbital = {}
        self.sum_correction_percentual = self._get_sum_correction_percentual()

        for symbol, orbitals in self.correction_indexes.items():
            for orbital in orbitals:
                cut = self._find_best_correction(symbol, orbital)
                cuts_per_atom_orbital[(symbol, orbital)] = cut

        gap = self._get_result_gap(is_indirect=self.indirect)
        return (cuts_per_atom_orbital, gap)

    def _find_best_correction(self, symbol: str, orbital: str) -> float:
        """
        Make calculation folder for correction of specific atom and orbital and find the best cut
            Args:
                symbol (str): Atom symbol
                orbital (str): Orbital type (s,p,d,f)
            Returns:
                gap_and_cut (tuple): Tuple containing
                the optimum cut and the gap generated
                by the correction.
        """
        folder_name = f"mkpotcar_{symbol.lower()}_{orbital.lower()}"
        path = os.path.join(self.root_folder, folder_name)
        if os.path.exists(path):
            shutil.rmtree(path)
        os.mkdir(path)

        number_equal_neighbors = self.divide_character[(symbol.capitalize(),
                                                        orbital.lower())]
        value = (100 / (1 + number_equal_neighbors)) * (
            self.band_projection[orbital][symbol] /
            self.sum_correction_percentual)
        logger.info(f"percentual of half electron is {round(value)}")

        cut = self._find_cut(symbol=symbol, base_path=path, orbital=orbital)

        potential_filename = get_output_filenames(
            "QE",
            self.input_files[0],
            atom=symbol
        )["potential"]

        # Copy corrected potential to the main folder
        find_cut_path = os.path.join(path, "find_cut")
        src = os.path.join(find_cut_path, "cut_{:.2f}".format(cut), os.path.basename(potential_filename))
        dest = os.path.dirname(os.path.abspath(self.input_files[0]))
        shutil.copy2(src, dest)

        return cut

    def _get_result_gap(self, is_indirect: bool) -> float:
        """
        Run a final ab initio calculation with all corrected UPF files
        and return the band gap.
        """
        calculation_folder = os.path.dirname(os.path.abspath(self.input_files[0]))

        self.runner.run(calculation_folder)

        filenames = get_output_filenames('QE', self.input_files[0])

        # All factory calls pass input_files to resolve prefix/outdir
        eigenvalues = self.software_factory.get_eigenvalues(
            base_path=calculation_folder,
            filename=filenames["eigenvalues"])
        fermi_energy = self.software_factory.get_fermi_energy(
            base_path=calculation_folder,
            filename=filenames["fermi_energy"])
        atoms_map = self.software_factory.get_atoms_map(
            base_path=calculation_folder,
            filename=filenames["atoms_map"])
        num_bands = self.software_factory.get_number_of_bands(
            base_path=calculation_folder,
            filename=filenames["number_of_bands"])
        band_projection_file = self.software_factory.get_band_projection_class(
            base_path=calculation_folder,
            filename=filenames["band_projection"])

        band_structure = BandStructure(eigenvalues, fermi_energy, atoms_map,
                                       num_bands, band_projection_file)
        gap_report = band_structure.band_gap(is_indirect)
        return gap_report["gap"]

    def _find_cut(self, symbol: str, base_path: str, orbital: str) -> float:
        """
        Find the cut which gives the maximum gap using Nelder-Mead.
        """

        folder = os.path.join(base_path, "find_cut")
        if os.path.exists(folder):
            shutil.rmtree(folder)
        os.mkdir(folder)

        function_args = {
            "base_path":                  base_path,
            "software_factory":           self.software_factory,
            "runner":                     self.runner,
            "symbol":                     symbol,
            "default_potential_filename": self.potential_filename,
            "atom_potential":             self.atom_potential,
            "potfiles_folder":            self.potential_folder,
            "amplitude":                  self.amplitude,
            "atoms":                      self.atoms,
            "software_files":             self.input_files,
            "is_conduction":              self.is_conduction,
            "indirect":                   self.indirect,
            "orbital":                    orbital,
            "exchange_correlation_type":  self.exchange_correlation_type,
            "max_iterations":             self.max_iterations,
            "calculation_code":           self.calculation_code,
            "ld1_command":                self.commands["ld1_command"],
            "virtual_v2_command":         self.commands["virtual_v2_command"]
        }
        cut_initial_guess = self.cut_initial_guess[(symbol.capitalize(),
                                                    orbital.lower())]
        result = minimize(
            find_negative_band_gap_qe,
            x0=cut_initial_guess,
            args=(function_args),
            method="Nelder-Mead",
            options={'xatol': self.tolerance}
        )
        if not result.success:
            logger.error("Optimization failed")
            raise Exception("Optimization failed.")

        return result.x[0]
    
    def _get_sum_correction_percentual(self) -> float:
        """
        Sum of the percentuals of orbitals to be corrected.
        """
        total_sum = 0
        for symbol, orbitals in self.correction_indexes.items():
            for orbital in orbitals:
                total_sum += self.band_projection[orbital][symbol]

        if total_sum == 0:
            logger.error(
                "No orbital selected for correction. Check your threshold.")
            raise ValueError(
                "No orbital selected for correction. Check your threshold.")

        return total_sum
