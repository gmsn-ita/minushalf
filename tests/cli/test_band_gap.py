"""
Test band gap command
"""
from click.testing import CliRunner
from minushalf.commands.band_gap import band_gap


def test_band_gap_gan_3d(file_path):
    """
    Test the result of band-gap for
    GaN 3d
    """
    base_path = file_path("/gan-3d/")
    result_path = file_path("/gan-3d/result_band_gap.txt")
    runner = CliRunner()
    result = runner.invoke(band_gap, ['-b', base_path])

    with open(result_path) as file:
        assert file.read() == result.output


def test_band_gap_bn_2d(file_path):
    """
    Test the result of band-gap for
    BN 2d
    """
    base_path = file_path("/bn-2d/")
    result_path = file_path("/bn-2d/result_band_gap.txt")
    runner = CliRunner()
    result = runner.invoke(band_gap, ['-b', base_path])

    with open(result_path) as file:
        assert file.read() == result.output


def test_band_gap_sic_2d(file_path):
    """
    Test the result of band-gap for
    SiC 2d
    """
    base_path = file_path("/sic-2d/")
    result_path = file_path("/sic-2d/result_band_gap.txt")
    runner = CliRunner()
    result = runner.invoke(band_gap, ['-b', base_path])

    with open(result_path) as file:
        assert file.read() == result.output


def test_band_gap_gec_2d(file_path):
    """
    Test the result of band-gap for
    GeC 2d
    """
    base_path = file_path("/gec-2d/")
    result_path = file_path("/gec-2d/result_band_gap.txt")
    runner = CliRunner()
    result = runner.invoke(band_gap, ['-b', base_path])

    with open(result_path) as file:
        assert file.read() == result.output


def test_band_gap_aln_2d(file_path):
    """
    Test the result of band-gap for
    AlN 2d
    """
    base_path = file_path("/aln-2d/")
    result_path = file_path("/aln-2d/result_band_gap.txt")
    runner = CliRunner()
    result = runner.invoke(band_gap, ['-b', base_path])

    with open(result_path) as file:
        assert file.read() == result.output
