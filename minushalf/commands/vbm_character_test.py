"""
Test vbm character command
"""
from click.testing import CliRunner
from minushalf.commands.vbm_character import vbm_character


def test_vbm_character_gan_3d(file_path):
    """
    Test the result of vbm-character for
    GaN 3d
    """
    base_path = file_path("/gan-3d/")
    result_path = file_path("/gan-3d/result_vbm_character.txt")
    runner = CliRunner()
    result = runner.invoke(vbm_character, ['-b', base_path])

    with open(result_path) as file:
        assert file.read() == result.output


def test_vbm_character_bn_2d(file_path):
    """
    Test the result of vbm-character for
    BN 2d
    """
    base_path = file_path("/bn-2d/")
    result_path = file_path("/bn-2d/result_vbm_character.txt")
    runner = CliRunner()
    result = runner.invoke(vbm_character, ['-b', base_path])

    with open(result_path) as file:
        assert file.read() == result.output


def test_vbm_character_sic_2d(file_path):
    """
    Test the result of vbm-character for
    SiC 2d
    """
    base_path = file_path("/sic-2d/")
    result_path = file_path("/sic-2d/result_vbm_character.txt")
    runner = CliRunner()
    result = runner.invoke(vbm_character, ['-b', base_path])

    with open(result_path) as file:
        assert file.read() == result.output


def test_vbm_character_gec_2d(file_path):
    """
    Test the result of vbm-character for
    GeC 2d
    """
    base_path = file_path("/gec-2d/")
    result_path = file_path("/gec-2d/result_vbm_character.txt")
    runner = CliRunner()
    result = runner.invoke(vbm_character, ['-b', base_path])

    with open(result_path) as file:
        assert file.read() == result.output


def test_vbm_character_aln_2d(file_path):
    """
    Test the result of vbm-character for
    AlN 2d
    """
    base_path = file_path("/aln-2d/")
    result_path = file_path("/aln-2d/result_vbm_character.txt")
    runner = CliRunner()
    result = runner.invoke(vbm_character, ['-b', base_path])

    with open(result_path) as file:
        assert file.read() == result.output
