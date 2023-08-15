"""
Trimming function
"""
import numpy as np
from minushalf.utils.constants import Constants


def trimming_function(
    radius: np.array,
    ion_potential: np.array,
    atom_potential: np.array,
    cut: float,
    amplitude: float,
) -> np.array:
    r"""
        Function that generate the potential for fractional occupation. The potential
        is cuted by a a function theta(r) to avoid divergence in calculations.
        The function of potential is defined as follows:

        .. math::
            V_{1/2} = (V_{atom}- V_{ion})\cdot \theta (r)

        where theta is:


        .. math::
            \theta (r)=A\cdot (1-(\frac{r}{CUT})^{n})^{3},r\leq CUT


        .. math::
            \theta (r) = 0,  r > CUT

        ---

            Args:
                cut (float): cutting parameter to cancel the potential

                amplitude (float): multiplicative factor of the potential function

                radius (np.array): rays in which the potential was calculated

                ion_potential (np.array): Atom pseudopotential with fractional occupation

                atom_potential (np.array): Atom pseudopotential with all electrons

            Returns:
                potential of fractional electron occupation at the exact level to be corrected


        """
    const = Constants()

    potential = np.array(
        4 * const.pi_constant * const.rydberg *
        np.power(const.bohr_radius, 3) *
        np.power(1 - np.power(radius / cut, const.trimming_exponent), 3) *
        (ion_potential - atom_potential) * amplitude)

    potential[radius >= cut] = 0

    return potential
