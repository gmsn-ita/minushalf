"""
Test band character command
"""
from click.testing import CliRunner
from minushalf.commands.band_character import band_character


def test_band_character_gan_3d(file_path):
    """
    Test the result of band-character for
    GaN 3d
    """
    base_path = file_path("/gan-3d/")
    result_path = file_path("/gan-3d/result_band_character.txt")
    runner = CliRunner()
    result = runner.invoke(band_character, ['4', '5', '-b', base_path])

    with open(result_path) as file:
        assert file.read() == result.output


def test_band_character_bn_2d(file_path):
    """
    Test the result of band-character for
    BN 2d
    """
    base_path = file_path("/bn-2d/")
    result_path = file_path("/bn-2d/result_band_character.txt")
    runner = CliRunner()
    result = runner.invoke(band_character, ['4', '5', '-b', base_path])

    with open(result_path) as file:
        assert file.read() == result.output


def test_band_character_sic_2d(file_path):
    """
    Test the result of band-character for
    SiC 2d
    """
    base_path = file_path("/sic-2d/")
    result_path = file_path("/sic-2d/result_band_character.txt")
    runner = CliRunner()
    result = runner.invoke(band_character, ['4', '5', '-b', base_path])

    with open(result_path) as file:
        assert file.read() == result.output


def test_band_character_gec_2d(file_path):
    """
    Test the result of band-character for
    GeC 2d
    """
    base_path = file_path("/gec-2d/")
    result_path = file_path("/gec-2d/result_band_character.txt")
    runner = CliRunner()
    result = runner.invoke(band_character, ['4', '5', '-b', base_path])

    with open(result_path) as file:
        assert file.read() == result.output


def test_band_character_aln_2d(file_path):
    """
    Test the result of band-character for
    AlN 2d
    """
    base_path = file_path("/aln-2d/")
    result_path = file_path("/aln-2d/result_band_character.txt")
    runner = CliRunner()
    result = runner.invoke(band_character, ['4', '5', '-b', base_path])

    with open(result_path) as file:
        assert file.read() == result.output
