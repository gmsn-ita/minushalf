"""
Test read minushalf yaml handler
"""
from minushalf.utils import minushalf_yaml
from minushalf.handlers import ReadMinushalf, read_minushalf_yaml
from minushalf.utils import MinushalfYaml

def test_default_parameters():
    """
    Test reading minushalf with the default parameters
    """
    read_minushalf = ReadMinushalf()

    request = read_minushalf.handle({})

    minushalf_yaml = request['minushalf_yaml']

    assert isinstance(minushalf_yaml,MinushalfYaml)
    assert minushalf_yaml.software == 'VASP'
    assert minushalf_yaml.correction['correction_code'] == 'v'
