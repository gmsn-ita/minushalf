"""
Fourier Transform
"""
import numpy as np
from minushalf.data import Constants


def correct_potential_fourier_transform(
    coefficient: np.array,
    k: np.array,
    rays: np.array,
    occupation_potential: np.array,
    cut: float,
) -> np.array:
    r"""
    The pseudopotential is given in terms of the radial distance, and is only defined for r >= 0,
    as expected. Since it is only evaluated inside an integral from 0 to infinity, it does not
    matter what values it assumes for r < 0. A natural choice is to define the function to be
    zero for negative values, but a more convenient choice is to choose
    v(-r)=-v(r) and n(-r)=-n(r), since purely real and odd functions have purely imaginary Fourier transforms.
    Let v' and n' be the odd extensions of the potential and the number density, respectively.

    .. math::
        E_{v} = \int_{0}^{\infty}v(r)n(r)dr = \int_{0}^{\infty}v'(r)n'(r)dr =
        \frac{1}{2}\cdot\int_{-\infty}^{\infty}v'(r)n'(r)dr = -\frac{1}{2}\cdot\int_{-\infty}^{\infty}V(k)N(k)dk


    On the third equalitty, we used the fact that the product of two odd
    functions is even, and in the last step we have applied Parseval's theorem,
    considering that the Fourier transforms are purely imaginary. Even though
    the function may not pass through the origin, we can still make an odd extension, by making it discontinuous.

    The data stored on POTCAR corresponds to the Fourier transform of the odd extension of v. It can be approximated
    by the summation on the right, where the prefactors were ommited.

    .. math::
        V(k) = i\cdot \sqrt{\frac{2}{\pi}}\cdot\int_{0}^{\infty}v(r)sin(b\cdot k\cdot r)dr
        \Rightarrow V(k)\sim \sum^{N_{r}}_{i=1}
        \frac{(v[i]\cdot sin(b\cdot k\cdot r[i])+v[i-1]\cdot sin(b\cdot k\cdot r[i-1]))}{2\cdot (r[i]-r[i-1])}

    Computes the opposite of the imaginary part of the j-th fourier transform coefficient
    through numerical integrationIndex zero stands for the r=DeltaR, and the function is
    assumed to be zero at the origin. Thus, the first trapezium of the numerical integration is
    degenerated to a triangulum, and its area must be calculated as so.

        Args:
            coefficient (np.array): Fourier transform of the potential for the atom in its ground state

            k (np.array): The wave vector in reciprocal space

            rays(np.array): List of rays on which pseudopotential calculations were made

            occupation_potential (np.arraygit): Potential of fractional electron occupation at the exact level to be corrected

            cut(float): Cutting parameter to cancel the potential

        Returns:
            Fourier transform of the potential for the state with fractional occupation of the crystal
    """
    const = Constants()
    k = np.reshape(k, (len(k), 1))

    if not k[0][0]:
        k[0][0] = 10**(-12)

    try:
        filter_rays = rays[np.where(rays < cut)]
    except ValueError as cut_error:
        raise ValueError(
            "the cut is smaller than all the rays passed") from cut_error

    radius = filter_rays[1:].astype(float)
    lazy_radius = filter_rays[:-1].astype(float)
    potential = occupation_potential[1:len(filter_rays)].astype(float)
    lazy_potential = occupation_potential[:len(filter_rays) - 1].astype(float)

    partial_fourier_sum = (
        (radius - lazy_radius) *
        (potential * np.sin(const.bohr_radius * k * radius) +
         lazy_potential * np.sin(const.bohr_radius * k * lazy_radius)) / 2)

    ## Sum all rows and transpose

    initial_term = (occupation_potential[0] *
                    np.sin(const.bohr_radius * k * rays[0]) / 2 * rays[0])

    ## Reshape vectors to give the correct output format
    initial_term = np.reshape(initial_term, (len(initial_term)))
    k = np.reshape(k, (len(k)))

    return (coefficient + (initial_term + partial_fourier_sum.sum(axis=1)) /
            (const.bohr_radius * k))
