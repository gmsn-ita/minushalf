"""
Function to check if a file exists
"""
import os


def check_vasprun_exists(func):
    """
    Function decrator to check if a file exists
    """
    def func_wrapper(self,
                     filename: str = "vasprun.xml",
                     base_path: str = None):
        path = filename
        if base_path:
            path = os.path.join(base_path, filename)
        if not os.path.exists(path):
            raise ValueError("File {} does not exist".format(filename))
        return func(self, filename, base_path)

    return func_wrapper


def check_procar_exists(func):
    """
    Function decrator to check if a file exists
    """
    def func_wrapper(self, filename: str = "PROCAR", base_path: str = None):
        path = filename
        if base_path:
            path = os.path.join(base_path, filename)
        if not os.path.exists(path):
            raise ValueError("File {} does not exist".format(filename))
        return func(self, filename, base_path)

    return func_wrapper


def check_eigenval_exists(func):
    """
    Function decrator to check if a file exists
    """
    def func_wrapper(self, filename: str = "EIGENVAL", base_path: str = None):
        path = filename
        if base_path:
            path = os.path.join(base_path, filename)
        if not os.path.exists(path):
            raise ValueError("File {} does not exist".format(filename))
        return func(self, filename, base_path)

    return func_wrapper


def check_potcar_exists(func):
    """
    Function decrator to check if a file exists
    """
    def func_wrapper(self, filename: str = "POTCAR", base_path: str = None):
        path = filename
        if base_path:
            path = os.path.join(base_path, filename)
        if not os.path.exists(path):
            raise ValueError("File {} does not exist".format(filename))
        return func(self, filename, base_path)

    return func_wrapper


def check_outcar_exists(func):
    """
    Function decrator to check if a file exists
    """
    def func_wrapper(self, *args, **kwargs):
        kwargs["filename"] = "OUTCAR"
        path = kwargs["filename"]
        if "base_path" in kwargs:
            path = os.path.join(kwargs["base_path"], kwargs["filename"])
        if not os.path.exists(path):
            raise ValueError("File {} does not exist".format(
                kwargs["filename"]))
        return func(self, *args, **kwargs)

    return func_wrapper
