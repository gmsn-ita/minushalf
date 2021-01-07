"""
Test vasp runner class
"""
import subprocess
from unittest import mock
from minushalf.softwares.vasp import VaspRunner


def test_with_number_of_cores_1():
    """
    Test vasprunner class with 1 core
    and path to vasp with default
    """
    with mock.patch("subprocess.run"):
        runner = VaspRunner()
        runner.run()
        command = ["mpirun", "-np", "1", "vasp"]
        subprocess.run.assert_called_once_with(command, check=True, cwd=".")


def test_with_number_of_cores_4():
    """
    Test vasprunner class with 4 core
    and path to vasp with '../vasp'
    """
    with mock.patch("subprocess.run"):
        runner = VaspRunner(path_to_vasp="../vasp", number_of_cores=4)
        runner.run()
        command = ["mpirun", "-np", "4", "../vasp"]
        subprocess.run.assert_called_once_with(command, check=True, cwd=".")


def test_with_different_path():
    """
    Test vasprunner class with a different execution path
    """
    with mock.patch("subprocess.run"):
        runner = VaspRunner(path_to_vasp="../vasp", number_of_cores=4)
        runner.run(cwd="./ss")
        command = ["mpirun", "-np", "4", "../vasp"]
        subprocess.run.assert_called_once_with(command, check=True, cwd="./ss")
