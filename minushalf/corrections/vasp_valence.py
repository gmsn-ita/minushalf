"""
Implements the algorithm that automates the process
of simple valence correction and optimizes the necessary parameters.
"""
import os
from subprocess import Popen, PIPE
import pandas as pd
from minushalf.atomic import InputFile


class ValenceCorrection():
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

    def execute(self):
        """
        Execute valence correction algorithm
        """
        atom_index = self._get_correction_atom_index()
        self._generate_atom_pseudopotential(atom_index)
        self._generate_occupation_potential(atom_index)
        #cut = self._find_cut(atom_index)

    def _find_cut(self):
        """
        Find the cut which gives the maximum gap
        """

    def _generate_occupation_potential(self, atom_index):
        """
        Generate the pseudo potential for the occupation
        of minus one half electron in the most relevant orbital
        on the VBM projection
        """
        atom_symbol = atom_index[0]
        secondary_quantum_number = atom_index[1]
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
            raise Exception("Call to occupation failed")

    def _generate_atom_pseudopotential(self, atom_index: tuple) -> None:
        """
        Make a dir with the atoms name,generate
        the input file and run the atomic program
        """
        atom_symbol = atom_index[0]
        folder_path = os.path.join(self.root_folder, atom_symbol)
        if os.path.exists(folder_path):
            os.remove(folder_path)
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
