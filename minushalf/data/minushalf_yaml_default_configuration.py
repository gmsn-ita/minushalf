"""
It helds the parameters and the
default values to be used in the
minushalf.yaml file
"""
from enum import Enum, unique


@unique
class CorrectionDefaultParams(Enum):
    """
    Default params for correction field
    """

    correction_code = "v"

    def __str__(self):
        return str(self.name)

    @staticmethod
    def to_list():
        """
        Return a list with default parameters
        """
        return list(map(lambda element: element.value,
                        CorrectionDefaultParams))

    @staticmethod
    def to_dict():
        """
        Returns a dictionary with the default parameters
        """
        values = map(lambda element: element.value, CorrectionDefaultParams)
        keys = map(lambda element: element.__str__(), CorrectionDefaultParams)
        return dict(zip(keys, values))


@unique
class VaspDefaultParams(Enum):
    """
    Default params for vasp
    """

    number_of_cores = 4
    path = "vasp"

    def __str__(self):
        return str(self.name)

    @staticmethod
    def to_list():
        """
        Return a list with default parameters
        """
        return list(map(lambda element: element.value, VaspDefaultParams))

    @staticmethod
    def to_dict():
        """
        Returns a dictionary with the default parameters
        """
        values = map(lambda element: element.value, VaspDefaultParams)
        keys = map(lambda element: element.__str__(), VaspDefaultParams)
        return dict(zip(keys, values))


@unique
class AtomicProgramDefaultParams(Enum):
    """
    Default params for atomic_program
    """

    exchange_correlation_code = "pb"
    calculation_code = "ae"
    max_iterations = 100

    def __str__(self):
        return str(self.name)

    @staticmethod
    def to_list():
        """
        Return a list with default parameters
        """
        return list(
            map(lambda element: element.value, AtomicProgramDefaultParams))

    @staticmethod
    def to_dict():
        """
        Returns a dictionary with the default parameters
        """
        values = map(lambda element: element.value, AtomicProgramDefaultParams)
        keys = map(lambda element: element.__str__(),
                   AtomicProgramDefaultParams)
        return dict(zip(keys, values))
