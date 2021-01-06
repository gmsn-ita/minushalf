"""
Trimming function
"""
import numpy as np
from minushalf.data import Constants


def trimming_function(radius: float, ion_potential: float, atom_potential,
                      cut: float, amplitude: float) -> float:
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
