"""
Trimming function
"""
import numpy as np
from minushalf.utils import Constants


def trimming_function(radius: float, ion_potential: float, atom_potential,
                      cut: float, amplitude: float) -> float:
    """
        Function that generate the potential for fractional occupation. The potential
        is cuted by a a function 0(r) to avoid divergence in calculations.
        The function of potential is defined as follows:

            V1/2 = (Vatom - Vion)*0(r), where 0(r) is;

                0(r) = A(1- (r/CUT)^n)^3 , for r <= CUT
                0(r) = 0, for r > CUT


            Args:
                cut (float): cutting parameter to cancel the potential

                amplitude (float): multiplicative factor of the potential function

                radius (float): radius in which the potential was calculated

                ion_potential (float): Atom pseudopotential with fractional occupation

                atom_potential (float): Atom pseudopotential with all electrons

            Returns:
                potential of fractional electron occupation at the exact level to be corrected


        """
    const = Constants()

    if radius >= cut:
        return 0

    return (4 * const.pi_constant * const.rydberg *
            np.power(const.bohr_radius, 3) *
            np.power(1 - np.power(radius / cut, const.trimming_exponent), 3) *
            (ion_potential - atom_potential) * amplitude)
