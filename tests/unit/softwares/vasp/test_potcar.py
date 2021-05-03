"""
Test Potcar class
"""
import numpy as np
from minushalf.softwares.vasp import Potcar


def test_potcar_ag(file_path):
    """
    Test the potcar parser for Silver
    """
    path = file_path('/Ag/POTCAR')
    potcar = Potcar(path)

    assert potcar.get_name() == "POTCAR"
    assert potcar.psctr_parameters[-1].strip() == "local part"
    assert potcar.psctr_parameters[3].strip() == "PAW_PBE Ag 34Apr4331"
    assert np.isclose(potcar.get_maximum_module_wave_vector(),
                      20.0022020200042)
    assert np.isclose(potcar.get_potential_fourier_transform()[0], 4.44444444)
    assert np.isclose(potcar.get_potential_fourier_transform()[-1], 3.33333333)
    assert potcar.last_lines[0].strip() == "gradient corrections used for XC"
    assert potcar.last_lines[-1].strip() == "End of Dataset"


def test_potcar_stringlist_ag(file_path):
    """
    Test the potcar function
    to_stringlist for Silver
    """
    path = file_path('/Ag/POTCAR')
    potcar = Potcar(path)

    with open(path, "r") as file:
        potcar_generated_lines = potcar.to_stringlist()
        for index, line in enumerate(file):
            assert line.strip() == potcar_generated_lines[index].strip()


def test_potcar_c(file_path):
    """
    Test the potcar parser for Carbon
    """
    path = file_path('/C/POTCAR')
    potcar = Potcar(path)

    assert potcar.get_name() == "POTCAR"
    assert potcar.psctr_parameters[-1].strip() == "local part"
    assert potcar.psctr_parameters[3].strip() == "PAW_PBE C 41Apr1441"
    assert np.isclose(potcar.get_maximum_module_wave_vector(),
                      232.232322323323)
    assert np.isclose(potcar.get_potential_fourier_transform()[0], 4.44444444)
    assert np.isclose(potcar.get_potential_fourier_transform()[-1], 2.22222222)
    assert potcar.last_lines[0].strip() == "gradient corrections used for XC"
    assert potcar.last_lines[-1].strip() == "End of Dataset"


def test_potcar_stringlist_c(file_path):
    """
    Test the potcar function
    to_stringlist for Carbon
    """
    path = file_path('/C/POTCAR')
    potcar = Potcar(path)

    with open(path, "r") as file:
        potcar_generated_lines = potcar.to_stringlist()
        for index, line in enumerate(file):
            assert line.strip() == potcar_generated_lines[index].strip()


def test_potcar_er(file_path):
    """
    Test the potcar parser for Erbium
    """
    path = file_path('/Er/POTCAR')
    potcar = Potcar(path)

    assert potcar.get_name() == "POTCAR"
    assert potcar.psctr_parameters[-1].strip() == "local part"
    assert potcar.psctr_parameters[3].strip() == "PAW_PBE Er 31Sep4331"
    assert np.isclose(potcar.get_maximum_module_wave_vector(),
                      23.0323023432223)
    assert np.isclose(potcar.get_potential_fourier_transform()[0], 2.22222222)
    assert np.isclose(potcar.get_potential_fourier_transform()[-1], 3.33333333)
    assert potcar.last_lines[0].strip() == "gradient corrections used for XC"
    assert potcar.last_lines[-1].strip() == "End of Dataset"


def test_potcar_stringlist_er(file_path):
    """
    Test the potcar function
    to_stringlist for Erbium
    """
    path = file_path('/Er/POTCAR')
    potcar = Potcar(path)

    with open(path, "r") as file:
        potcar_generated_lines = potcar.to_stringlist()
        for index, line in enumerate(file):
            assert line.strip() == potcar_generated_lines[index].strip()


def test_potcar_f(file_path):
    """
    Test the potcar parser for Fluorine
    """
    path = file_path('/F/POTCAR')
    potcar = Potcar(path)

    assert potcar.get_name() == "POTCAR"
    assert potcar.psctr_parameters[-1].strip() == "local part"
    assert potcar.psctr_parameters[3].strip() == "PAW_PBE F 11Apr1111"
    assert np.isclose(potcar.get_maximum_module_wave_vector(),
                      144.114114424213)
    assert np.isclose(potcar.get_potential_fourier_transform()[0], 3.34333333)
    assert np.isclose(potcar.get_potential_fourier_transform()[-1], 1.11111111)
    assert potcar.last_lines[0].strip() == "gradient corrections used for XC"
    assert potcar.last_lines[-1].strip() == "End of Dataset"


def test_potcar_stringlist_f(file_path):
    """
    Test the potcar function
    to_stringlist for Fluorine
    """
    path = file_path('/F/POTCAR')
    potcar = Potcar(path)

    with open(path, "r") as file:
        potcar_generated_lines = potcar.to_stringlist()
        for index, line in enumerate(file):
            assert line.strip() == potcar_generated_lines[index].strip()


def test_potcar_fe(file_path):
    """
    Test the potcar parser for Iron
    """
    path = file_path('/Fe/POTCAR')
    potcar = Potcar(path)

    assert potcar.get_name() == "POTCAR"
    assert potcar.psctr_parameters[-1].strip() == "local part"
    assert potcar.psctr_parameters[3].strip() == "PAW_PBE Fe 21Sep2222"
    assert np.isclose(potcar.get_maximum_module_wave_vector(),
                      22.2222222222142)
    assert np.isclose(potcar.get_potential_fourier_transform()[0], 4.34344444)
    assert np.isclose(potcar.get_potential_fourier_transform()[-1], 4.41442111)
    assert potcar.last_lines[0].strip() == "gradient corrections used for XC"
    assert potcar.last_lines[-1].strip() == "End of Dataset"


def test_potcar_stringlist_fe(file_path):
    """
    Test the potcar function
    to_stringlist for Iron
    """
    path = file_path('/Fe/POTCAR')
    potcar = Potcar(path)

    with open(path, "r") as file:
        potcar_generated_lines = potcar.to_stringlist()
        for index, line in enumerate(file):
            assert line.strip() == potcar_generated_lines[index].strip()


def test_potcar_ga(file_path):
    """
    Test the potcar parser for Gallium
    """
    path = file_path('/Ga/POTCAR')
    potcar = Potcar(path)

    assert potcar.get_name() == "POTCAR"
    assert potcar.psctr_parameters[-1].strip() == "local part"
    assert potcar.psctr_parameters[3].strip() == "PAW_PBE Ga 33Apr3333"
    assert np.isclose(potcar.get_maximum_module_wave_vector(),
                      22.3000320300002)
    assert np.isclose(potcar.get_potential_fourier_transform()[0], 4.44444444)
    assert np.isclose(potcar.get_potential_fourier_transform()[-1], 2.22222224)
    assert potcar.last_lines[0].strip() == "gradient corrections used for XC"
    assert potcar.last_lines[-1].strip() == "End of Dataset"


def test_potcar_stringlist_ga(file_path):
    """
    Test the potcar function
    to_stringlist for Gallium
    """
    path = file_path('/Ga/POTCAR')
    potcar = Potcar(path)

    with open(path, "r") as file:
        potcar_generated_lines = potcar.to_stringlist()
        for index, line in enumerate(file):
            assert line.strip() == potcar_generated_lines[index].strip()


def test_potcar_h(file_path):
    """
    Test the potcar parser for Hydrogen
    """
    path = file_path('/H/POTCAR')
    potcar = Potcar(path)

    assert potcar.get_name() == "POTCAR"
    assert potcar.psctr_parameters[-1].strip() == "local part"
    assert potcar.psctr_parameters[3].strip() == "PAW_PBE H 33Jun3333"
    assert np.isclose(potcar.get_maximum_module_wave_vector(),
                      121.120000122110)
    assert np.isclose(potcar.get_potential_fourier_transform()[0], 3.33333333)
    assert np.isclose(potcar.get_potential_fourier_transform()[-1], 4.44444444)
    assert potcar.last_lines[0].strip() == "gradient corrections used for XC"
    assert potcar.last_lines[-1].strip() == "End of Dataset"


def test_potcar_stringlist_h(file_path):
    """
    Test the potcar function
    to_stringlist for Hydrogen
    """
    path = file_path('/H/POTCAR')
    potcar = Potcar(path)

    with open(path, "r") as file:
        potcar_generated_lines = potcar.to_stringlist()
        for index, line in enumerate(file):
            assert line.strip() == potcar_generated_lines[index].strip()


def test_potcar_he(file_path):
    """
    Test the potcar parser for Helium
    """
    path = file_path('/He/POTCAR')
    potcar = Potcar(path)

    assert potcar.get_name() == "POTCAR"
    assert potcar.psctr_parameters[-1].strip() == "local part"
    assert potcar.psctr_parameters[3].strip() == "PAW_PBE He 33Jan4334"
    assert np.isclose(potcar.get_maximum_module_wave_vector(),
                      022.224440022024)
    assert np.isclose(potcar.get_potential_fourier_transform()[0], 3.33333333)
    assert np.isclose(potcar.get_potential_fourier_transform()[-1], 4.44444444)
    assert potcar.last_lines[0].strip() == "gradient corrections used for XC"
    assert potcar.last_lines[-1].strip() == "End of Dataset"


def test_potcar_stringlist_he(file_path):
    """
    Test the potcar function
    to_stringlist for Helium
    """
    path = file_path('/He/POTCAR')
    potcar = Potcar(path)

    with open(path, "r") as file:
        potcar_generated_lines = potcar.to_stringlist()
        for index, line in enumerate(file):
            assert line.strip() == potcar_generated_lines[index].strip()


def test_potcar_mg(file_path):
    """
    Test the potcar parser for Magnesium
    """
    path = file_path('/Mg/POTCAR')
    potcar = Potcar(path)

    assert potcar.get_name() == "POTCAR"
    assert potcar.psctr_parameters[-1].strip() == "local part"
    assert potcar.psctr_parameters[3].strip() == "PAW_PBE Mg 11Apr1112"
    assert np.isclose(potcar.get_maximum_module_wave_vector(),
                      43.3004443404304)
    assert np.isclose(potcar.get_potential_fourier_transform()[0], 1.11111311)
    assert np.isclose(potcar.get_potential_fourier_transform()[-1], 3.33333332)
    assert potcar.last_lines[0].strip() == "gradient corrections used for XC"
    assert potcar.last_lines[-1].strip() == "End of Dataset"


def test_potcar_stringlist_mg(file_path):
    """
    Test the potcar function
    to_stringlist for Magnesium
    """
    path = file_path('/Mg/POTCAR')
    potcar = Potcar(path)

    with open(path, "r") as file:
        potcar_generated_lines = potcar.to_stringlist()
        for index, line in enumerate(file):
            assert line.strip() == potcar_generated_lines[index].strip()


def test_potcar_ir(file_path):
    """
    Test the potcar parser for Iridium
    """
    path = file_path('/Ir/POTCAR')
    potcar = Potcar(path)

    assert potcar.get_name() == "POTCAR"
    assert potcar.psctr_parameters[-1].strip() == "local part"
    assert potcar.psctr_parameters[3].strip() == "PAW_PBE Ir 13Sep1111"
    assert np.isclose(potcar.get_maximum_module_wave_vector(),
                      32.3233333333332)
    assert np.isclose(potcar.get_potential_fourier_transform()[0], 4.44444444)
    assert np.isclose(potcar.get_potential_fourier_transform()[-1], 4.41114111)
    assert potcar.last_lines[0].strip() == "gradient corrections used for XC"
    assert potcar.last_lines[-1].strip() == "End of Dataset"


def test_potcar_stringlist_ir(file_path):
    """
    Test the potcar function
    to_stringlist for Iridium
    """
    path = file_path('/Ir/POTCAR')
    potcar = Potcar(path)

    with open(path, "r") as file:
        potcar_generated_lines = potcar.to_stringlist()
        for index, line in enumerate(file):
            assert line.strip() == potcar_generated_lines[index].strip()


def test_potcar_mo(file_path):
    """
    Test the potcar parser for Molybdenum
    """
    path = file_path('/Mo/POTCAR')
    potcar = Potcar(path)

    assert potcar.get_name() == "POTCAR"
    assert potcar.psctr_parameters[-1].strip() == "local part"
    assert potcar.psctr_parameters[3].strip() == "PAW_PBE Mo 33Apr1331"
    assert np.isclose(potcar.get_maximum_module_wave_vector(),
                      22.2220200222202)
    assert np.isclose(potcar.get_potential_fourier_transform()[0], 2.22222222)
    assert np.isclose(potcar.get_potential_fourier_transform()[-1], 1.11111111)
    assert potcar.last_lines[0].strip() == "gradient corrections used for XC"
    assert potcar.last_lines[-1].strip() == "End of Dataset"


def test_potcar_stringlist_mo(file_path):
    """
    Test the potcar function
    to_stringlist for Molybdenum
    """
    path = file_path('/Mo/POTCAR')
    potcar = Potcar(path)

    with open(path, "r") as file:
        potcar_generated_lines = potcar.to_stringlist()
        for index, line in enumerate(file):
            assert line.strip() == potcar_generated_lines[index].strip()


def test_potcar_n(file_path):
    """
    Test the potcar parser for Nitrogen
    """
    path = file_path('/N/POTCAR')
    potcar = Potcar(path)

    assert potcar.get_name() == "POTCAR"
    assert potcar.psctr_parameters[-1].strip() == "local part"
    assert potcar.psctr_parameters[3].strip() == "PAW_PBE N 44Apr4444"
    assert np.isclose(potcar.get_maximum_module_wave_vector(),
                      434.434441343344)
    assert np.isclose(potcar.get_potential_fourier_transform()[0], 4.33444344)
    assert np.isclose(potcar.get_potential_fourier_transform()[-1], 2.22342222)
    assert potcar.last_lines[0].strip() == "gradient corrections used for XC"
    assert potcar.last_lines[-1].strip() == "End of Dataset"


def test_potcar_stringlist_n(file_path):
    """
    Test the potcar function
    to_stringlist for Nitrogen
    """
    path = file_path('/N/POTCAR')
    potcar = Potcar(path)

    with open(path, "r") as file:
        potcar_generated_lines = potcar.to_stringlist()
        for index, line in enumerate(file):
            assert line.strip() == potcar_generated_lines[index].strip()


def test_potcar_ne(file_path):
    """
    Test the potcar parser for Neon
    """
    path = file_path('/Ne/POTCAR')
    potcar = Potcar(path)

    assert potcar.get_name() == "POTCAR"
    assert potcar.psctr_parameters[-1].strip() == "local part"
    assert potcar.psctr_parameters[3].strip() == "PAW_PBE Ne 43Jan2443"
    assert np.isclose(potcar.get_maximum_module_wave_vector(),
                      220.424220000201)
    assert np.isclose(potcar.get_potential_fourier_transform()[0], 3.33333333)
    assert np.isclose(potcar.get_potential_fourier_transform()[-1], 4.44444444)
    assert potcar.last_lines[0].strip() == "gradient corrections used for XC"
    assert potcar.last_lines[-1].strip() == "End of Dataset"


def test_potcar_stringlist_ne(file_path):
    """
    Test the potcar function
    to_stringlist for Neon
    """
    path = file_path('/Ne/POTCAR')
    potcar = Potcar(path)

    with open(path, "r") as file:
        potcar_generated_lines = potcar.to_stringlist()
        for index, line in enumerate(file):
            assert line.strip() == potcar_generated_lines[index].strip()


def test_potcar_ru(file_path):
    """
    Test the potcar parser for Ruthenium
    """
    path = file_path('/Ru/POTCAR')
    potcar = Potcar(path)

    assert potcar.get_name() == "POTCAR"
    assert potcar.psctr_parameters[-1].strip() == "local part"
    assert potcar.psctr_parameters[3].strip() == "PAW_PBE Ru 31Feb0334"
    assert np.isclose(potcar.get_maximum_module_wave_vector(),
                      14.4341134441414)
    assert np.isclose(potcar.get_potential_fourier_transform()[0], 3.33333333)
    assert np.isclose(potcar.get_potential_fourier_transform()[-1], 4.44444444)
    assert potcar.last_lines[0].strip() == "gradient corrections used for XC"
    assert potcar.last_lines[-1].strip() == "End of Dataset"


def test_potcar_stringlist_ru(file_path):
    """
    Test the potcar function
    to_stringlist for Ruthenium
    """
    path = file_path('/Ru/POTCAR')
    potcar = Potcar(path)

    with open(path, "r") as file:
        potcar_generated_lines = potcar.to_stringlist()
        for index, line in enumerate(file):
            assert line.strip() == potcar_generated_lines[index].strip()


def test_potcar_s(file_path):
    """
    Test the potcar parser for Sulfur
    """
    path = file_path('/S/POTCAR')
    potcar = Potcar(path)

    assert potcar.get_name() == "POTCAR"
    assert potcar.psctr_parameters[-1].strip() == "local part"
    assert potcar.psctr_parameters[3].strip() == "PAW_PBE S 43Sep3444"
    assert np.isclose(potcar.get_maximum_module_wave_vector(),
                      23.2322222332323)
    assert np.isclose(potcar.get_potential_fourier_transform()[0], 4.44444444)
    assert np.isclose(potcar.get_potential_fourier_transform()[-1], 4.44444444)
    assert potcar.last_lines[0].strip() == "gradient corrections used for XC"
    assert potcar.last_lines[-1].strip() == "End of Dataset"


def test_potcar_stringlist_s(file_path):
    """
    Test the potcar function
    to_stringlist for Sulfur
    """
    path = file_path('/S/POTCAR')
    potcar = Potcar(path)

    with open(path, "r") as file:
        potcar_generated_lines = potcar.to_stringlist()
        for index, line in enumerate(file):
            assert line.strip() == potcar_generated_lines[index].strip()


def test_potcar_sn(file_path):
    """
    Test the potcar parser for Tin
    """
    path = file_path('/Sn/POTCAR')
    potcar = Potcar(path)

    assert potcar.get_name() == "POTCAR"
    assert potcar.psctr_parameters[-1].strip() == "local part"
    assert potcar.psctr_parameters[3].strip() == "PAW_PBE Sn 13Apr3113"
    assert np.isclose(potcar.get_maximum_module_wave_vector(),
                      44.4444234442442)
    assert np.isclose(potcar.get_potential_fourier_transform()[0], 1.11111111)
    assert np.isclose(potcar.get_potential_fourier_transform()[-1], 1.41441111)
    assert potcar.last_lines[0].strip() == "gradient corrections used for XC"
    assert potcar.last_lines[-1].strip() == "End of Dataset"


def test_potcar_stringlist_sn(file_path):
    """
    Test the potcar function
    to_stringlist for Tin
    """
    path = file_path('/Sn/POTCAR')
    potcar = Potcar(path)

    with open(path, "r") as file:
        potcar_generated_lines = potcar.to_stringlist()
        for index, line in enumerate(file):
            assert line.strip() == potcar_generated_lines[index].strip()


def test_potcar_tc(file_path):
    """
    Test the potcar parser for Technetium
    """
    path = file_path('/Tc/POTCAR')
    potcar = Potcar(path)

    assert potcar.get_name() == "POTCAR"
    assert potcar.psctr_parameters[-1].strip() == "local part"
    assert potcar.psctr_parameters[3].strip() == "PAW_PBE Tc 10Feb0114"
    assert np.isclose(potcar.get_maximum_module_wave_vector(),
                      11.1111111111111)
    assert np.isclose(potcar.get_potential_fourier_transform()[0], 2.22222222)
    assert np.isclose(potcar.get_potential_fourier_transform()[-1], 1.41121414)
    assert potcar.last_lines[0].strip() == "gradient corrections used for XC"
    assert potcar.last_lines[-1].strip() == "End of Dataset"


def test_potcar_stringlist_tc(file_path):
    """
    Test the potcar function
    to_stringlist for Technetium
    """
    path = file_path('/Tc/POTCAR')
    potcar = Potcar(path)

    with open(path, "r") as file:
        potcar_generated_lines = potcar.to_stringlist()
        for index, line in enumerate(file):
            assert line.strip() == potcar_generated_lines[index].strip()


def test_potcar_v(file_path):
    """
    Test the potcar parser for Vanadium
    """
    path = file_path('/V/POTCAR')
    potcar = Potcar(path)

    assert potcar.get_name() == "POTCAR"
    assert potcar.psctr_parameters[-1].strip() == "local part"
    assert potcar.psctr_parameters[3].strip() == "PAW_PBE V 11Apr3113"
    assert np.isclose(potcar.get_maximum_module_wave_vector(),
                      33.3331333111131)
    assert np.isclose(potcar.get_potential_fourier_transform()[0], 1.11111111)
    assert np.isclose(potcar.get_potential_fourier_transform()[-1], 2.22444222)
    assert potcar.last_lines[0].strip() == "gradient corrections used for XC"
    assert potcar.last_lines[-1].strip() == "End of Dataset"


def test_potcar_stringlist_v(file_path):
    """
    Test the potcar function
    to_stringlist for Vanadium
    """
    path = file_path('/V/POTCAR')
    potcar = Potcar(path)

    with open(path, "r") as file:
        potcar_generated_lines = potcar.to_stringlist()
        for index, line in enumerate(file):
            assert line.strip() == potcar_generated_lines[index].strip()


def test_potcar_w(file_path):
    """
    Test the potcar parser for Tungsten
    """
    path = file_path('/W/POTCAR')
    potcar = Potcar(path)

    assert potcar.get_name() == "POTCAR"
    assert potcar.psctr_parameters[-1].strip() == "local part"
    assert potcar.psctr_parameters[3].strip() == "PAW_PBE W 41Apr4444"
    assert np.isclose(potcar.get_maximum_module_wave_vector(),
                      42.4444422222444)
    assert np.isclose(potcar.get_potential_fourier_transform()[0], 4.44444444)
    assert np.isclose(potcar.get_potential_fourier_transform()[-1], 4.44444444)
    assert potcar.last_lines[0].strip() == "gradient corrections used for XC"
    assert potcar.last_lines[-1].strip() == "End of Dataset"


def test_potcar_stringlist_w(file_path):
    """
    Test the potcar function
    to_stringlist for Tungsten
    """
    path = file_path('/W/POTCAR')
    potcar = Potcar(path)

    with open(path, "r") as file:
        potcar_generated_lines = potcar.to_stringlist()
        for index, line in enumerate(file):
            assert line.strip() == potcar_generated_lines[index].strip()


def test_potcar_xe(file_path):
    """
    Test the potcar parser for Xenon
    """
    path = file_path('/Xe/POTCAR')
    potcar = Potcar(path)

    assert potcar.get_name() == "POTCAR"
    assert potcar.psctr_parameters[-1].strip() == "local part"
    assert potcar.psctr_parameters[3].strip() == "PAW_PBE Xe 21Sep2222"
    assert np.isclose(potcar.get_maximum_module_wave_vector(), 3.3220323033312)
    assert np.isclose(potcar.get_potential_fourier_transform()[0], 1.11111111)
    assert np.isclose(potcar.get_potential_fourier_transform()[-1], 4.44424424)
    assert potcar.last_lines[0].strip() == "gradient corrections used for XC"
    assert potcar.last_lines[-1].strip() == "End of Dataset"


def test_potcar_stringlist_xe(file_path):
    """
    Test the potcar function
    to_stringlist for Xenon
    """
    path = file_path('/Xe/POTCAR')
    potcar = Potcar(path)

    with open(path, "r") as file:
        potcar_generated_lines = potcar.to_stringlist()
        for index, line in enumerate(file):
            assert line.strip() == potcar_generated_lines[index].strip()


def test_potcar_lda_si(file_path):
    """
    Test the potcar using LDA for silicium
    """
    path = file_path('/Si/POTCAR_LDA')
    potcar = Potcar(path)

    assert potcar.get_name() == "POTCAR"
    assert potcar.psctr_parameters[-1].strip() == "local part"
    assert potcar.psctr_parameters[3].strip() == "PAW Si 04Apr4222"
    assert np.isclose(potcar.get_maximum_module_wave_vector(),
                      31.2131311111111)
    assert np.isclose(potcar.get_potential_fourier_transform()[0], 2.22222222)
    assert np.isclose(potcar.get_potential_fourier_transform()[-1], 3.33333333)
    assert potcar.last_lines[0].strip() == "core charge-density (partial)"
    assert potcar.last_lines[-1].strip() == "End of Dataset"


def test_potcar_stringlist_si_lda(file_path):
    """
    Test the potcar function
    to_stringlist for Silicium. The
    potcar uses LDA approximation.
    """
    path = file_path('/Si/POTCAR_LDA')
    potcar = Potcar(path)

    with open(path, "r") as file:
        potcar_generated_lines = potcar.to_stringlist()
        for index, line in enumerate(file):
            assert line.strip() == potcar_generated_lines[index].strip()
