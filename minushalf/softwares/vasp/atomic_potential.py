"""
Correct crystal potential to
fractional occupations
"""
import copy
import numpy as np
from minushalf.atomic import Vtotal, InputFile
from minushalf.utils import (trimming_function,
                             correct_potential_fourier_transform)


class AtomicPotential():
    """
    Correct atomic potential fourier tranform
    for fractional occupations in valence
    or conduction bands
    """
    def __init__(self, vtotal: Vtotal, vtotal_occupied: Vtotal,
                 input_file: InputFile, potcar: any) -> None:
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
        self.input_file = input_file
        self.potcar = potcar

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
        wave_vectors = np.arange((len(self.potcar.fourier_coefficients))(
            self.potcar.k_max / len(self.potcar.fourier_coefficients)))

        correct_potential = np.vectorize(correct_potential_fourier_transform)
        potential = correct_potential(
            self.potcar.fourier_coefficients,
            wave_vectors,
            self.vtotal.radius,
            occupation_potential,
            cut,
        )

        return potential

    def correct_file(self, potential: list, cut: float,
                     amplitude: float) -> None:
        """
        Create the potential file corrected

            Args:
                potential (list): List of corrected potentials fourier transform

                cut (float): Cutting parameter to cancel the potential

                amplitude (float): Multiplicative factor of the potential function
        """

        filename = ""
        if np.isclose(abs(amplitude), 1.0):
            filename = "POTCARcut{:.2f}".format(cut)
        else:
            filename = "POTCARcut{:.2f}A{:.1f}".format(cut, amplitude)

        copy_potcar = copy.deepcopy(self.potcar)
        copy_potcar.potential = potential
        copy_potcar.to_file(filename)
