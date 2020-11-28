"""
Test file for projected_wave_function module
"""
from minushalf.projected_wave_function import ProjectedWaveFunction
from minushalf.projected_wave_function.vasp import Vasp


def test_factory_vasp(file_path):
    """
    Test if object created by the
    factory is an instance of Vasp class.
    """
    path_to_vasprun = file_path("vasprun_GaN.xml")

    factory = ProjectedWaveFunction()
    vasp = factory.vasp(path_to_vasprun.__str__())

    assert isinstance(vasp, Vasp)
