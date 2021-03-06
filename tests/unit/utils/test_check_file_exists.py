"""
Test to check if file exists
"""
import pytest
from minushalf.utils import (
    check_potcar_exists,
    check_eigenval_exists,
    check_vasprun_exists,
    check_procar_exists,
    check_outcar_exists,
)


def test_potcar_exists(file_path):
    """
    Test function if potcar exists
    """
    potcar_path = file_path("/Ag/")
    func = lambda _, filename, base_path=None: f"{filename}{base_path}"
    func_wrapper = check_potcar_exists(func)
    func_wrapper({}, base_path=potcar_path)


@pytest.mark.xfail
def test_potcar_dont_exists(file_path):
    """
    Test function if potcar exists
    """
    potcar_path = file_path("/Co/")
    func = lambda _, filename, base_path=None: f"{filename}{base_path}"
    func_wrapper = check_potcar_exists(func)
    func_wrapper({}, base_path=potcar_path)


def test_vasprun_exists(file_path):
    """
    Test function if vasprun exists
    """
    potcar_path = file_path("/gec-2d/")
    func = lambda _, filename, base_path=None: f"{filename}{base_path}"
    func_wrapper = check_vasprun_exists(func)
    func_wrapper({}, base_path=potcar_path)


@pytest.mark.xfail
def test_vasprun_dont_exists(file_path):
    """
    Test function if vasprun exists
    """
    potcar_path = file_path("/Co/")
    func = lambda _, filename, base_path=None: f"{filename}{base_path}"
    func_wrapper = check_vasprun_exists(func)
    func_wrapper({}, base_path=potcar_path)


def test_procar_exists(file_path):
    """
    Test function if procar exists
    """
    potcar_path = file_path("/gec-2d/")
    func = lambda _, filename, base_path=None: f"{filename}{base_path}"
    func_wrapper = check_procar_exists(func)
    func_wrapper({}, base_path=potcar_path)


@pytest.mark.xfail
def test_procar_dont_exists(file_path):
    """
    Test function if procar exists
    """
    potcar_path = file_path("/Co/")
    func = lambda _, filename, base_path=None: f"{filename}{base_path}"
    func_wrapper = check_procar_exists(func)
    func_wrapper({}, base_path=potcar_path)


def test_eigenval_exists(file_path):
    """
    Test function if eigenval exists
    """
    potcar_path = file_path("/gec-2d/")
    func = lambda _, filename, base_path=None: f"{filename}{base_path}"
    func_wrapper = check_eigenval_exists(func)
    func_wrapper({}, base_path=potcar_path)


@pytest.mark.xfail
def test_eigenval_dont_exists(file_path):
    """
    Test function if eigenval exists
    """
    potcar_path = file_path("/Co/")
    func = lambda _, filename, base_path=None: f"{filename}{base_path}"
    func_wrapper = check_eigenval_exists(func)
    func_wrapper({}, base_path=potcar_path)


def test_outcar_exists(file_path):
    """
    Test function if outcar file exists
    """
    outcar_path = file_path("/sic-2d/")
    func = lambda _, ion_index, filename, base_path=None: f"{filename}{base_path}"
    func_wrapper = check_outcar_exists(func)
    func_wrapper({}, '', base_path=outcar_path)


@pytest.mark.xfail
def test_outcar_dont_exists(file_path):
    """
    Test function if outcar file exists
    """
    outcar_path = file_path("/Co/")
    func = lambda _, ion_index, filename, base_path=None: f"{filename}{base_path}"
    func_wrapper = check_outcar_exists(func)
    func_wrapper({}, '', base_path=outcar_path)
