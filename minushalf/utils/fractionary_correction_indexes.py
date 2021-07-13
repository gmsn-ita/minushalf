"""
Get atoms to be corrected for
simple valence correction and simple
condunction correction
"""
import itertools
import pandas as pd
from collections import defaultdict


def get_fractionary_correction_indexes(band_projection: pd.DataFrame,
                                        treshold: int = 5) -> dict:
    """
        Get dataframe index of the orbitals which contributes more
        than 5 percent to (VBM|CBM)

            Returns:

                correction_indexes (dict):A dict wherw the keys are the atoms
                                          symbols and the value is a list with the orbitals type to be
                                          corrected.
                                          Ex:
                                          {
                                          'Ga': ['p','s'],
                                          'N' : ['d','f'],
                                          }

    """
    columns_name = band_projection.columns.tolist()
    rows_name = band_projection.index.tolist()

    correction_indexes = defaultdict(list)
    for row, column in itertools.product(rows_name, columns_name):
        if band_projection[column][row] >= treshold:
            correction_indexes[row].append(column)
    return correction_indexes
