"""
Makes ternary search and find the
maximum of a function unimodal
"""
import numpy as np


def ternary_search(minimum: float, maximum: float,
                   cost_function: any) -> float:
    """
    Given a unimodal cost function, it realizes ternary
    search and return the maximum value of the function
        Args:
            minimum (float): minimum value of the search
            maximum (float): maximum value of the search
            cost_function (function): function to be evaluated
        Returns:
            x (float ): value at which the function reaches its maximum
    """
    while not np.isclose(maximum, minimum):
        mid_left = minimum + (maximum - minimum) / 3
        mid_right = maximum - (maximum - minimum) / 3
        result_mid_left = cost_function(mid_left)
        result_mid_right = cost_function(mid_right)
        if result_mid_left < result_mid_right:
            minimum = mid_left
        else:
            maximum = mid_right

    return minimum
