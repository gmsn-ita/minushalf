"""
Test software configurations input class
"""
from minushalf.io import SoftwareConfigurations


def test_default_parameters():
    """
    Test the default parameters 
    """
    software_configurations = SoftwareConfigurations()

    assert software_configurations.software_name == "VASP"
    assert software_configurations.command == ["mpirun", "vasp"]


def test_override_parameters():
    """
    Test the class with other parameters
    """
    software_configurations = SoftwareConfigurations(
        software_name="VASP", command=["mpirun", "-np", "4"])

    assert software_configurations.software_name == "VASP"
    assert software_configurations.command == ["mpirun", "-np", "4"]


def test_to_list():
    """
    Test method to_list
    """
    software_configurations = SoftwareConfigurations()

    parameters_list = software_configurations.to_list()

    assert parameters_list == [["mpirun", "vasp"]]


def test_to_dict():
    """
    Test method to_list
    """
    software_configurations = SoftwareConfigurations()

    parameters_dict = software_configurations.to_dict()

    assert parameters_dict == {
        "command": ["mpirun", "vasp"],
    }
