"""
Test band character command
"""
from click.testing import CliRunner
from minushalf.commands.band_character import band_character


def test_band_character_gan_3d_VASP(file_path):
    """
    Test the result of band-character for
    GaN 3d in VASP format.
    """
    base_path = file_path("/gan-3d/vasp/")
    result_path = file_path("/gan-3d/vasp/result_band_character.txt")
    runner = CliRunner()
    result = runner.invoke(band_character, ['4', '5', '-b', base_path])

    with open(result_path) as file:
        assert file.read() == result.output


def test_band_character_bn_2d_VASP(file_path):
    """
    Test the result of band-character for
    BN 2d in VASP format.
    """
    base_path = file_path("/bn-2d/")
    result_path = file_path("/bn-2d/result_band_character.txt")
    runner = CliRunner()
    result = runner.invoke(band_character, ['4', '5', '-b', base_path])

    with open(result_path) as file:
        assert file.read() == result.output


def test_band_character_sic_2d_VASP(file_path):
    """
    Test the result of band-character for
    SiC 2d in VASP format.
    """
    base_path = file_path("/sic-2d/")
    result_path = file_path("/sic-2d/result_band_character.txt")
    runner = CliRunner()
    result = runner.invoke(band_character, ['4', '5', '-b', base_path])

    with open(result_path) as file:
        assert file.read() == result.output


def test_band_character_gec_2d_VASP(file_path):
    """
    Test the result of band-character for
    GeC 2d in VASP format.
    """
    base_path = file_path("/gec-2d/")
    result_path = file_path("/gec-2d/result_band_character.txt")
    runner = CliRunner()
    result = runner.invoke(band_character, ['4', '5', '-b', base_path])

    with open(result_path) as file:
        assert file.read() == result.output


def test_band_character_aln_2d_VASP(file_path):
    """
    Test the result of band-character for
    AlN 2d in VASP format.
    """
    base_path = file_path("/aln-2d/")
    result_path = file_path("/aln-2d/result_band_character.txt")
    runner = CliRunner()
    result = runner.invoke(band_character, ['4', '5', '-b', base_path])

    with open(result_path) as file:
        assert file.read() == result.output

def test_band_character_gan_3d_QE(file_path):
    """
    Test the result of band-character for
    GaN 3d in QE format.
    """
    base_path = file_path("/gan-3d/qe/")
    input_path = file_path("/gan-3d/qe/proj.in")
    result_path = file_path("/gan-3d/qe/result_band_character_qe.txt")
    runner = CliRunner()
    result = runner.invoke(band_character, ['6', '7', '-b', base_path, '-s', 'QE', '-n', input_path])

    with open(result_path) as file:
        assert file.read() == result.output