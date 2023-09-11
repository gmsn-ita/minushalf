"""
Fourier Transform
"""
import numpy as np
from minushalf.utils.constants import Constants


def correct_potential_fourier_transform(
    coefficient: np.array,
    k: np.array,
    rays: np.array,
    occupation_potential: np.array,
    cut: float,
) -> np.array:
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

    # Sum all rows and transpose

    initial_term = (occupation_potential[0] *
                    np.sin(const.bohr_radius * k * rays[0]) / 2 * rays[0])

    # Reshape vectors to give the correct output format
    initial_term = np.reshape(initial_term, (len(initial_term)))
    k = np.reshape(k, (len(k)))

    return (coefficient + (initial_term + partial_fourier_sum.sum(axis=1)) /
            (const.bohr_radius * k))
