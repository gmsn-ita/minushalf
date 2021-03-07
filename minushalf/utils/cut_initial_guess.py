"""
Gives cut initial guess
"""
from minushalf.data import CutInitialGuessMethods


def _3d_crystals(distance: float) -> float:
    """
        Gives cut initial guess for 3D crystals

            Args:
                distance (float): Nearest neighbor distance.

            Returns:
                cut_guess (float): An initial guess to cut.
        """
    return 0.15 + 0.84 * distance


class CutInitialGuess:
    """
    Estimate cut inital guess from
    the nearest neighbor distance.
    """
    def __init__(self):
        """
        Describe options to guess
        """
        self.options = {
            CutInitialGuessMethods.three_dimensions.value: _3d_crystals,
        }

    def guess(self, distance: float, method: str) -> float:
        """
        Given the nearest neighbor distance and the method, it returns the
        initial guess.

            Args:
                method (str): method of gessing
                distance (float): nearest neighbor distance

            Returns:
                cut_guess (float): An initial guess to cut.
        """
        return self.options[method](distance)
