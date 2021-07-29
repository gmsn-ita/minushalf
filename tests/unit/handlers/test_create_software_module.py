"""
Test create software module handler
"""
from minushalf.handlers import CreateSoftwareModule
from minushalf.io import MinushalfYaml
from minushalf.softwares import Vasp


def test_if_object_was_created():
    """
    Test if the object was realy created
    """
    minushalf_yaml = MinushalfYaml.from_file()
    request = {"minushalf_yaml": minushalf_yaml}
    handler = CreateSoftwareModule()

    response = handler.action(request)

    assert isinstance(response["software_module"], Vasp)
