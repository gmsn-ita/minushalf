"""
Correct crystal potential to
fractional occupations
"""
import copy
import numpy as np
from minushalf.softwares.potential_file import PotentialFile
from minushalf.io.vtotal import Vtotal
from minushalf.utils.trimming_function import trimming_function
from minushalf.utils.correct_potential_fourier_transform import correct_potential_fourier_transform


class AtomicPotential():
    """
    Correct atomic potential fourier tranform
    for fractional occupations in valence
    or conduction bands
    """
    def __init__(
        self,
        vtotal: Vtotal,
        vtotal_occupied: Vtotal,
        potential_file: PotentialFile,
    ) -> None:
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
        trimmed_half_electron_potential = trimming_function(
            np.array(self.vtotal.radius, dtype=float),
            np.array(self.vtotal_occupied.down_potential, dtype=float),
            np.array(self.vtotal.down_potential, dtype=float),
            cut,
            amplitude,
        )
        return trimmed_half_electron_potential

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

        occupied_potential = self.occupy_potential(cut, amplitude)

        ## generate a sample of the absolute values of the wave vector
        size_potential_sampling = len(
            self.potential_file.get_potential_fourier_transform())
        maximum_absolute_value_wave_vector = self.potential_file.get_maximum_module_wave_vector(
        )
        absolute_wave_vectors = np.arange(size_potential_sampling) * (
            maximum_absolute_value_wave_vector / size_potential_sampling) + (maximum_absolute_value_wave_vector / size_potential_sampling)

        corrected_potential = correct_potential_fourier_transform(
            coefficient=np.array(
                self.potential_file.get_potential_fourier_transform()),
            k=absolute_wave_vectors,
            rays=np.array(self.vtotal.radius, dtype=float),
            occupation_potential=np.array(occupied_potential, dtype=float),
            cut=cut,
        )

        return corrected_potential

    def correct_file(
        self,
        potential: list,
        cut: float,
        amplitude: float,
        is_conduction: bool = False,
    ) -> None:
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
        if np.isclose(amplitude, 1.0):
            filename = "{}cut{:.2f}".format(self.potential_file.get_name(),
                                            cut)
        else:
            filename = "{}cut{:.2f}A{:.1f}".format(
                self.potential_file.get_name(), cut, amplitude)

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
