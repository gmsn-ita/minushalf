"""
Fourier Transform
"""
import numpy as np
from minushalf.utils import Constants


def correct_potential_fourier_transform(coefficient: float, k: float,
                                        rays: np.array,
                                        occupation_potential: np.array,
                                        cut: float) -> float:
    r"""
    The pseudopotential is given in terms of the radial distance, and is only defined for r >= 0,
    as expected. Since it is only evaluated inside an integral from 0 to infinity, it does not
    matter what values it assumes for r < 0. A natural choice is to define the function to be
    zero for negative values, but a more convenient choice is to choose
    v(-r)=-v(r) and n(-r)=-n(r), since purely real and odd functions have purely imaginary Fourier transforms.
    Let v' and n' be the odd extensions of the potential and the number density, respectively.

         /\ inf          /\ inf                  /\ inf                   /\ inf
         |               |                       |                        |
    Ev = |  v(r)n(r)dr = |  v'(r)n'(r)dr = (1/2) |  v'(r)n'(r)dr = -(1/2) | V(k)N(k)dk
         |               |                       |                        |
        \/  0           \/  0                   \/ -inf                  \/ -inf
    On the third equalitty, we used the fact that the product of two odd
    functions is even, and in the last step we have applied Parseval's theorem,
    considering that the Fourier transforms are purely imaginary. Even though
    the function may not pass through the origin, we can still make an odd extension, by making it discontinuous.

    The data stored on POTCAR corresponds to the Fourier transform of the odd extension of v. It can be approximated
    by the summation on the right, where the prefactors were ommited.
                                                    ______
                        /\ inf                      \     | Nr
                        |                            \
    V(k) = i*sqrt(2/pi) |  v(r)sin(bkr)dr  => V(k) ~  >     (v[i]sin(bkr[i]) + v[i-1]sin(bkr[i-1]))/2 *(r[i]-r[i-1])
                        |                            /
                       \/ 0                         /_____| 1

    Computes the opposite of the imaginary part of the j-th fourier transform coefficient
    through numerical integrationIndex zero stands for the r=DeltaR, and the function is
    assumed to be zero at the origin. Thus, the first trapezium of the numerical integration is
    degenerated to a triangulum, and its area must be calculated as so.

        Args:
            coefficient (float): Fourier transform of the potential for the atom in its ground state

            k (float): The wave vector in reciprocal space

            rays(list): List of rays on which pseudopotential calculations were made

            occupation_potential (list): Potential of fractional electron occupation at the exact level to be corrected

            cut(float): Cutting parameter to cancel the potential

        Returns:
            Fourier transform of the potential for the state with fractional occupation of the crystal
    """
    const = Constants()

    if not k:
        k = 10**(-12)

    try:
        filter_rays = rays[np.where(rays < cut)]
    except ValueError as cut_error:
        raise ValueError(
            "the cut is smaller than all the rays passed") from cut_error

    partial_fourier_sum = lambda radius, lazy_radius, potential, lazy_potential: (
        (potential * np.sin(const.bohr_radius * k * radius) + lazy_potential *
         np.sin(const.bohr_radius * k * lazy_radius)) / 2 *
        (radius - lazy_radius))

    get_fourier_components = np.vectorize(partial_fourier_sum)

    fourier_components = get_fourier_components(
        filter_rays[1:],
        filter_rays[:-1],
        occupation_potential[1:len(filter_rays)],
        occupation_potential[:len(filter_rays) - 1],
    )

    initial_term = (occupation_potential[0] *
                    np.sin(const.bohr_radius * k * rays[0]) / 2 * rays[0])

    return (coefficient + (initial_term + fourier_components.sum()) /
            (const.bohr_radius * k))
