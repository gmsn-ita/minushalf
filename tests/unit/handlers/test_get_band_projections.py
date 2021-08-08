"""
Test get projections handler
"""
from unittest.mock import patch
from minushalf.io import MinushalfYaml
from minushalf.softwares import Vasp
from minushalf.handlers import GetBandProjections
from minushalf.utils import BandStructure


def test_get_valence_projections(file_path):
    """
    Get projections of GaN 3d
    """
    ## Variables
    minushalf_yaml = MinushalfYaml.from_file()
    minushalf_yaml.correction.correction_code = "vc"
    software_module = Vasp()
    request = {
        "minushalf_yaml": minushalf_yaml,
        "software_module": software_module,
    }
    band_structure = BandStructure.create(software_module,
                                          file_path("/gan-3d/"))

    ## Mock method in class BandStructure
    with patch('minushalf.utils.BandStructure.create') as mock:
        mock.return_value = band_structure
        handler = GetBandProjections()
        response = handler.action(request)

    ## Assertions
    assert response["projections"]["valence"]["d"]["Ga"] == 16
    assert response["projections"]["conduction"]["s"]["N"] == 56
    assert "projections" in response


def test_get_conduction_and_valence_projections(file_path):
    """
    Get projections of GaN 3d
    """
    ## Variables
    minushalf_yaml = MinushalfYaml.from_file()
    software_module = Vasp()
    request = {
        "minushalf_yaml": minushalf_yaml,
        "software_module": software_module,
    }
    band_structure = BandStructure.create(software_module,
                                          file_path("/gan-3d/"))

    ## Mock method in class BandStructure
    with patch('minushalf.utils.BandStructure.create') as mock:
        mock.return_value = band_structure
        handler = GetBandProjections()
        response = handler.action(request)

    ## Assertions
    assert response["projections"]["valence"]["d"]["Ga"] == 16
    assert response["projections"]["valence"]["p"]["N"] == 78
    assert "projections" in response


def test_get_projections_with_ovewrite(file_path):
    """
    Get overwrited projections of GaN 3d
    """
    ## Variables
    minushalf_yaml = MinushalfYaml.from_file(
        file_path("/minushalf_yaml/minushalf_filled_out.yaml"))
    software_module = Vasp()
    request = {
        "minushalf_yaml": minushalf_yaml,
        "software_module": software_module,
    }
    band_structure = BandStructure.create(software_module,
                                          file_path("/gan-3d/"))

    ## Mock method in class BandStructure
    with patch('minushalf.utils.BandStructure.create') as mock:
        mock.return_value = band_structure
        handler = GetBandProjections()
        response = handler.action(request)
    print(response["projections"])
    ## Assertions
    assert response["projections"]["valence"]["d"]["Ga"] == 98
    assert "projections" in response
    assert "valence" in response["projections"]
    assert ("conduction" in response["projections"]) == False
