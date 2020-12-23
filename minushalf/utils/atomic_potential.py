"""
Correct crystal potential to
fractional occupations
"""
import copy
import numpy as np
from minushalf.interfaces import PotentialFile
from .vtotal import Vtotal
from .trimming_function import trimming_function
from .correct_potential_fourier_transform import correct_potential_fourier_transform


class AtomicPotential():
    """
    Correct atomic potential fourier tranform
    for fractional occupations in valence
    or conduction bands
    """
    def __init__(self, vtotal: Vtotal, vtotal_occupied: Vtotal,
                 potential_file: PotentialFile) -> None:
        """
            Args:
                vtotal (Vtotal): Class of the atom's VTOTAL file
                with all electrons in the ground state

                vtotal_occupied (Vtotal):Class of the atom's VTOTAL file
                with fractional occupation

                input_file (InputFile): Class of the INP file used
                to generate the VTOTAL with fractional occupation

                potcar (Potcar): Class of the POTCAR file of the
                respective atom

        """
        self.vtotal = vtotal
        self.vtotal_occupied = vtotal_occupied
        self.potential_file = potential_file

    def occupy_potential(self, cut: float, amplitude) -> list:
        """
            Args:
                cut (float): Cutting parameter to cancel the potential

                amplitude (float): Multiplicative factor of the potential function

            Returns:
                A list that contains the potentials of fractional electron
                occupation at the exact level to be corrected.
        """
        trimming = np.vectorize(trimming_function)
        occupation_potential = trimming(
            self.vtotal.radius,
            self.vtotal_occupied.down_potential,
            self.vtotal.down_potential,
            cut,
            amplitude,
        )
        return occupation_potential

    def correct_potential(self,
                          cut: float,
                          amplitude: float,
                          is_conduction: bool = False) -> list:
        """
        Correct fourier transform of the potential (V(k)) present in POTCAR file.

            Args:
                cut (float): Cutting parameter to cancel the potential

                amplitude (float): Multiplicative factor of the potential function

                is_conduction (bool): Indicates whether the potential correction will be
                made in the valence or in the conduction

            Returns:
                List of corrected potentials fourier transform


        """
        if is_conduction:
            amplitude = abs(amplitude) * -1

        occupation_potential = self.occupy_potential(cut, amplitude)
        wave_vectors = np.arange(
            len(self.potential_file.get_potential_fourier_transform())) * (
                self.potential_file.get_maximum_module_wave_vector() /
                len(self.potential_file.get_potential_fourier_transform()))

        correct_potential = np.vectorize(
            correct_potential_fourier_transform,
            excluded=['rays', 'occupation_potential'],
        )

        potential = correct_potential(
            coefficient=self.potential_file.get_potential_fourier_transform(),
            k=wave_vectors,
            rays=np.array(self.vtotal.radius, dtype=object),
            occupation_potential=np.array(occupation_potential, dtype=object),
            cut=cut,
        )

        return potential

    def correct_file(self,
                     potential: list,
                     cut: float,
                     amplitude: float,
                     is_conduction: bool = False) -> None:
        """
        Create the potential file corrected

            Args:
                potential (list): List of corrected potentials fourier transform

                cut (float): Cutting parameter to cancel the potential

                amplitude (float): Multiplicative factor of the potential function
        """
        if is_conduction:
            amplitude = abs(amplitude) * -1

        filename = ""
        if np.isclose(abs(amplitude), 1.0):
            filename = "POTCARcut{:.2f}".format(cut)
        else:
            filename = "POTCARcut{:.2f}A{:.1f}".format(cut, amplitude)

        copy_potcar = copy.deepcopy(self.potential_file)
        copy_potcar.potential = potential
        copy_potcar.to_file(filename)

    def get_corrected_file_lines(
        self,
        potential: list,
    ) -> list:
        """
        Create the potential file corrected

            Args:
                potential (list): List of corrected potentials fourier transform
            Returns:
                potential_lines(list): A List of potcar lines
        """
        copy_potcar = copy.deepcopy(self.potential_file)
        copy_potcar.potential = potential
        return copy_potcar.to_stringlist()
