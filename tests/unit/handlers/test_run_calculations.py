"""
Test handler that runs first principle calculations
"""
from unittest.mock import MagicMock, Mock
from minushalf.handlers import RunCalculations
from minushalf.io import MinushalfYaml
from minushalf.softwares import Vasp


def test_run_calculations():
    """
    Check if the right functions are called
    """
    ## Creating mocks
    run_mock = Mock()
    runner_mock = Mock(run=run_mock)
    get_runner_mock = Mock(return_value=runner_mock)

    software_module = Vasp()
    software_module.get_runner = get_runner_mock

    request = {
        "minushalf_yaml": MinushalfYaml.from_file(),
        "software_module": software_module
    }

    handler = RunCalculations()
    handler.action(request)

    assert get_runner_mock.assert_called_once_with(['mpirun', 'vasp']) == None
    assert get_runner_mock.called
    assert run_mock.called
    assert run_mock.assert_called_once() == None
