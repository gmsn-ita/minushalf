"""
Test vbm character command
"""
from click.testing import CliRunner
from minushalf.commands import vbm_character


def test_vbm_character_gan_3d(file_path):
    """
    Test the result of vbm-character for
    GaN 3d
    """
    procar_filename = file_path("/gan-3d/PROCAR")
    eigenval_filename = file_path("/gan-3d/EIGENVAL")
    vasprun_filename = file_path("/gan-3d/vasprun.xml")
    result_path = file_path("/gan-3d/result_vbm_character.txt")
    runner = CliRunner()
    result = runner.invoke(vbm_character, [
        '-p', procar_filename, '-e', eigenval_filename, '-v', vasprun_filename
    ])

    with open(result_path) as file:
        assert file.read() == result.output
