"""
Test correct potfile command..
"""
import os
import numpy as np
from click.testing import CliRunner
from minushalf.commands import correct_potfile


def test_correct_potfile_ga_first(file_path):
    """
    Test correct_potfile command with Ga
    """
    runner = CliRunner()
    vtotal = file_path("/Ga/VTOTAL")
    vtotal_occ = file_path("/Ga/VTOTAL_OCC")
    inp_occ = file_path("/Ga/INP_OCC")
    potcar = file_path("/Ga/")

    cut = 3.56
    amplitude = 1.0
    with runner.isolated_filesystem():
        result = runner.invoke(correct_potfile, [
            "-C", "3.56", "-b", "{}".format(potcar), "-o",
            "{}".format(vtotal_occ), "-v", "{}".format(vtotal)
        ])
        assert result.exit_code == 0

        filename = ""
        if np.isclose(abs(amplitude), 1.0):
            filename = "POTCARcut{:.2f}".format(cut)
        else:
            filename = "POTCARcut{:.2f}A{:.1f}".format(cut, amplitude)

        assert os.path.exists(filename) == True


def test_correct_potfile_ga_second(file_path):
    """
    Test correct_potfile command with Ga
    """
    runner = CliRunner()
    vtotal = file_path("/Ga/VTOTAL")
    vtotal_occ = file_path("/Ga/VTOTAL_OCC")
    potcar = file_path("/Ga/")

    cut = 3.56
    amplitude = 1.2
    with runner.isolated_filesystem():
        result = runner.invoke(correct_potfile, [
            "-C", "3.56", "-b", "{}".format(potcar), "-o",
            "{}".format(vtotal_occ), "-v", "{}".format(vtotal), "-a",
            "{}".format(amplitude)
        ])
        assert result.exit_code == 0

        filename = ""
        if np.isclose(abs(amplitude), 1.0):
            filename = "POTCARcut{:.2f}".format(cut)
        else:
            filename = "POTCARcut{:.2f}A{:.1f}".format(cut, amplitude)

        assert os.path.exists(filename) == True


def test_correct_potfile_ga_third(file_path):
    """
    Test correct_potfile command with Ga
    """
    runner = CliRunner()
    vtotal = file_path("/Ga/VTOTAL")
    vtotal_occ = file_path("/Ga/VTOTAL_OCC")
    potcar = file_path("/Ga/")

    cuts = np.arange(3.5, 3.8, 0.1)
    amplitude = 1.2
    with runner.isolated_filesystem():
        result = runner.invoke(correct_potfile, [
            "-C", "3.5:0.1:3.8", "-b", "{}".format(potcar), "-o",
            "{}".format(vtotal_occ), "-v", "{}".format(vtotal), "-a",
            "{}".format(amplitude)
        ])
        assert result.exit_code == 0
        for cut in cuts:
            filename = ""
            if np.isclose(abs(amplitude), 1.0):
                filename = "POTCARcut{:.2f}".format(cut)
            else:
                filename = "POTCARcut{:.2f}A{:.1f}".format(cut, amplitude)

            assert os.path.exists(filename) == True


def test_correct_potfile_ga_fourth(file_path):
    """
    Test correct_potfile command with Ga
    """
    runner = CliRunner()
    vtotal = file_path("/Ga/VTOTAL")
    vtotal_occ = file_path("/Ga/VTOTAL_OCC")
    potcar = file_path("/Ga/")

    cuts = np.arange(3.5, 3.8, 0.1)
    amplitude = 1.0
    with runner.isolated_filesystem():
        result = runner.invoke(correct_potfile, [
            "-C", "3.5:0.1:3.8", "-b", "{}".format(potcar), "-o",
            "{}".format(vtotal_occ), "-v", "{}".format(vtotal)
        ])
        assert result.exit_code == 0
        for cut in cuts:
            filename = ""
            if np.isclose(abs(amplitude), 1.0):
                filename = "POTCARcut{:.2f}".format(cut)
            else:
                filename = "POTCARcut{:.2f}A{:.1f}".format(cut, amplitude)

            assert os.path.exists(filename) == True


def test_correct_potfile_ga_third_conduction_correction(file_path):
    """
    Test correct_potfile command with Ga with conduction correction
    """
    runner = CliRunner()
    vtotal = file_path("/Ga/VTOTAL")
    vtotal_occ = file_path("/Ga/VTOTAL_OCC")
    potcar = file_path("/Ga/")

    cuts = np.arange(3.5, 3.8, 0.1)
    amplitude = 1.2
    with runner.isolated_filesystem():
        result = runner.invoke(correct_potfile, [
            "-C", "3.5:0.1:3.8", "-b", "{}".format(potcar), "-o",
            "{}".format(vtotal_occ), "-v", "{}".format(vtotal), "-a",
            "{}".format(amplitude), "-c", "CONDUCTION"
        ])
        assert result.exit_code == 0
        for cut in cuts:
            filename = ""
            if np.isclose(abs(amplitude), 1.0):
                filename = "POTCARcut{:.2f}".format(cut)
            else:
                filename = "POTCARcut{:.2f}A-{:.1f}".format(cut, amplitude)

            assert os.path.exists(filename) == True
