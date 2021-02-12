"""
Parse cut string used in correct potential file
"""
import numpy as np


def parse_cut(cut: str) -> list:
    """
    Parse cut in a list of numbers.

    Args:
        cut (str): Cut energy to be used in the program, it can be
                   passed in two ways:

                    unique value : float or integer
                    range:  begin(float|integer):pass(float|integer):end(float|integer)
    Returns:

        cut_numbers (list): Permited values of cut.
    """
    if not cut:
        raise ValueError('A valeu of cut energy must be provided.')

    cuts = cut.split(':')

    try:
        cut_number = [float(element) for element in cuts]
    except ValueError as bad_input:
        raise ValueError("Invalid Input.") from bad_input

    if len(cut_number) == 1:
        return cut_number
    elif len(cut_number) != 3:
        raise ValueError()

    return np.arange(cut_number[0],
                     cut_number[2],
                     cut_number[1],
                     dtype=np.float64)
