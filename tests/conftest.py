"""
File for test configuration
"""
import os
import pytest


@pytest.fixture
def file_path():
    """
    Return the path to any fixture file.

    Args:
        filename (str): name of the file.
    """

    def _foo(filename: str):

        path_to_script = os.path.realpath(__file__)
        current_directory = os.path.split(path_to_script)[0]
        return os.path.join(current_directory, "fixtures/{}".format(filename)).__str__()

    return _foo
