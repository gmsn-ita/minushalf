"""
Transform informations about band_project generated
by cbm_character, vbm_character or band_character in
a normalized datafame Grouped by orbital types
"""
from collections import defaultdict
import pandas as pd
from minushalf.utils.orbital import Orbital, OrbitalType


def projection_to_df(projection: defaultdict(list)) -> pd.DataFrame:
    """
    Transform received dictionaries into information with a higher degree of readability.

        Args:
            projection (defaultdict(list)): A dictionary containing the projections
            of each atom per orbital
        Returns:
            prolection_df (pd.DataFrame): A dataframe containing the projections
            per orbital type and normalized between 0 - 100
    """
    orbitals = [orbital.__str__() for orbital in Orbital]
    number_columns = 0
    for value in list(projection.values()):
        number_columns = max(number_columns, len(value))

    if not number_columns:
        raise ValueError('The projection vector is empty')

    if number_columns > len(orbitals):
        raise ValueError(
            'the projection vector has more columns than necessary')

    projection_df = (pd.DataFrame.from_dict(projection,
                                            orient="index",
                                            columns=orbitals[:number_columns]))
    normalized_df = projection_df.div(projection_df.sum().sum())
    normalized_df = normalized_df.mul(100)

    orbital_type = [str(elem) for elem in OrbitalType]

    def join_orbitals(cols):
        return [
            cat for col in cols for cat in orbital_type if col.startswith(cat)
        ]

    normalized_df = normalized_df.groupby(join_orbitals(projection_df.columns),
                                          axis=1).sum().round()
    return normalized_df
