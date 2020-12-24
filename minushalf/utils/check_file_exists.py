"""
Function to check if a file exists
"""
import os


def check_file_exists(func):
    """
    Function decrator to check if a file exists
    """
    def func_wrapper(filename: str, base_path: str = None):
        if base_path:
            path = os.path.join(base_path, filename)
        if not os.path.exists(path):
            raise ValueError("File {} does not exist".format(filename))
        return func(filename, base_path)

    return func_wrapper
