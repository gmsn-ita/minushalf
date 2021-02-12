"""
Get atoms to be corrected for
simple valence correction and simple
condunction correction
"""
import itertools
import math
import pandas as pd
import numpy as np


def get_simple_corrections_indexes(band_projection: pd.DataFrame) -> dict:
    """
        Get dataframe index of the orbital which contributes more
        to (VBM|CBM)

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

    maximum_projection = -math.inf

    for row, column in itertools.product(rows_name, columns_name):
        projection = band_projection[column][row]
        maximum_projection = max(projection, maximum_projection)
        if np.isclose(projection, maximum_projection):
            max_row = row
            max_column = column

    correction_indexes = {max_row: [max_column]}
    return correction_indexes
