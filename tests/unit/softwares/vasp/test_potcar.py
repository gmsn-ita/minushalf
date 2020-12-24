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

    assert potcar.psctr_parameters[-1].strip() == "local part"
    assert potcar.psctr_parameters[0].strip() == "PAW_PBE Ag 02Apr2005"
    assert np.isclose(potcar.k_max, 75.5890395431569)
    assert np.isclose(potcar.fourier_coefficients[0], 0.77325279e2)
    assert np.isclose(potcar.fourier_coefficients[-1], 0.34828552e0)
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

    assert potcar.psctr_parameters[-1].strip() == "local part"
    assert potcar.psctr_parameters[0].strip() == "PAW_PBE C 08Apr2002"
    assert np.isclose(potcar.k_max, 124.721915246209)
    assert np.isclose(potcar.fourier_coefficients[0], 0.20536679e2)
    assert np.isclose(potcar.fourier_coefficients[-1], 0.46962343e-1)
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

    assert potcar.psctr_parameters[-1].strip() == "local part"
    assert potcar.psctr_parameters[0].strip() == "PAW_PBE Er 01Sep2006"
    assert np.isclose(potcar.k_max, 71.8095875659991)
    assert np.isclose(potcar.fourier_coefficients[0], 0.33161532e3)
    assert np.isclose(potcar.fourier_coefficients[-1], 0.77488015e0)
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

    assert potcar.psctr_parameters[-1].strip() == "local part"
    assert potcar.psctr_parameters[0].strip() == "PAW_PBE F 08Apr2002"
    assert np.isclose(potcar.k_max, 122.832189257630)
    assert np.isclose(potcar.fourier_coefficients[0], 0.25090202e2)
    assert np.isclose(potcar.fourier_coefficients[-1], 0.84227234e-1)
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

    assert potcar.psctr_parameters[-1].strip() == "local part"
    assert potcar.psctr_parameters[0].strip() == "PAW_PBE Fe 06Sep2000"
    assert np.isclose(potcar.k_max, 81.2582175088937)
    assert np.isclose(potcar.fourier_coefficients[0], 0.86392220e2)
    assert np.isclose(potcar.fourier_coefficients[-1], 0.21998356e0)
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

    assert potcar.psctr_parameters[-1].strip() == "local part"
    assert potcar.psctr_parameters[0].strip() == "PAW_PBE Ga 08Apr2002"
    assert np.isclose(potcar.k_max, 71.8095875659991)
    assert np.isclose(potcar.fourier_coefficients[0], 0.59460173e2)
    assert np.isclose(potcar.fourier_coefficients[-1], 0.10631365e0)
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

    assert potcar.psctr_parameters[-1].strip() == "local part"
    assert potcar.psctr_parameters[0].strip() == "PAW_PBE H 15Jun2001"
    assert np.isclose(potcar.k_max, 170.075338972103)
    assert np.isclose(potcar.fourier_coefficients[0], 0.24361675E1)
    assert np.isclose(potcar.fourier_coefficients[-1], 0.62759860e-2)
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

    assert potcar.psctr_parameters[-1].strip() == "local part"
    assert potcar.psctr_parameters[0].strip() == "PAW_PBE He 05Jan2001"
    assert np.isclose(potcar.k_max, 170.075338972103)
    assert np.isclose(potcar.fourier_coefficients[0], 0.48278121e1)
    assert np.isclose(potcar.fourier_coefficients[-1], 0.12539022e-1)
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

    assert potcar.psctr_parameters[-1].strip() == "local part"
    assert potcar.psctr_parameters[0].strip() == "PAW_PBE Mg 13Apr2007"
    assert np.isclose(potcar.k_max, 94.4862994289462)
    assert np.isclose(potcar.fourier_coefficients[0], 0.22502855e2)
    assert np.isclose(potcar.fourier_coefficients[-1], 0.40609408e-1)
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

    assert potcar.psctr_parameters[-1].strip() == "local part"
    assert potcar.psctr_parameters[0].strip() == "PAW_PBE Ir 06Sep2000"
    assert np.isclose(potcar.k_max, 71.8095875659991)
    assert np.isclose(potcar.fourier_coefficients[0], 0.12669216e3)
    assert np.isclose(potcar.fourier_coefficients[-1], 0.31623286e0)
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

    assert potcar.psctr_parameters[-1].strip() == "local part"
    assert potcar.psctr_parameters[0].strip() == "PAW_PBE Mo 08Apr2002"
    assert np.isclose(potcar.k_max, 68.0301355888413)
    assert np.isclose(potcar.fourier_coefficients[0], 0.11440305e3)
    assert np.isclose(potcar.fourier_coefficients[-1], 0.23531220e0)
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

    assert potcar.psctr_parameters[-1].strip() == "local part"
    assert potcar.psctr_parameters[0].strip() == "PAW_PBE N 08Apr2002"
    assert np.isclose(potcar.k_max, 124.721915246209)
    assert np.isclose(potcar.fourier_coefficients[0], 0.21549196e2)
    assert np.isclose(potcar.fourier_coefficients[-1], 0.58613255e-1)
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

    assert potcar.psctr_parameters[-1].strip() == "local part"
    assert potcar.psctr_parameters[0].strip() == "PAW_PBE Ne 05Jan2001"
    assert np.isclose(potcar.k_max, 109.604107337578)
    assert np.isclose(potcar.fourier_coefficients[0], 0.19440536e2)
    assert np.isclose(potcar.fourier_coefficients[-1], 0.12060531e0)
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

    assert potcar.psctr_parameters[-1].strip() == "local part"
    assert potcar.psctr_parameters[0].strip() == "PAW_PBE Ru 04Feb2005"
    assert np.isclose(potcar.k_max, 69.9198615774202)
    assert np.isclose(potcar.fourier_coefficients[0], 0.57352320e2)
    assert np.isclose(potcar.fourier_coefficients[-1], 0.29557935e0)
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

    assert potcar.psctr_parameters[-1].strip() == "local part"
    assert potcar.psctr_parameters[0].strip() == "PAW_PBE S 06Sep2000"
    assert np.isclose(potcar.k_max, 98.2657514061040)
    assert np.isclose(potcar.fourier_coefficients[0], 0.25157640e2)
    assert np.isclose(potcar.fourier_coefficients[-1], 0.11327915e0)
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

    assert potcar.psctr_parameters[-1].strip() == "local part"
    assert potcar.psctr_parameters[0].strip() == "PAW_PBE Sn 08Apr2002"
    assert np.isclose(potcar.k_max, 62.3609576231045)
    assert np.isclose(potcar.fourier_coefficients[0], 0.31273018e2)
    assert np.isclose(potcar.fourier_coefficients[-1], 0.18728585e0)
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

    assert potcar.psctr_parameters[-1].strip() == "local part"
    assert potcar.psctr_parameters[0].strip() == "PAW_PBE Tc 04Feb2005"
    assert np.isclose(potcar.k_max, 66.1404096002623)
    assert np.isclose(potcar.fourier_coefficients[0], 0.31304226e2)
    assert np.isclose(potcar.fourier_coefficients[-1], 0.29089674e0)
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

    assert potcar.psctr_parameters[-1].strip() == "local part"
    assert potcar.psctr_parameters[0].strip() == "PAW_PBE V 08Apr2002"
    assert np.isclose(potcar.k_max, 69.9198615774202)
    assert np.isclose(potcar.fourier_coefficients[0], 0.61993946e2)
    assert np.isclose(potcar.fourier_coefficients[-1], 0.18556083e0)
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

    assert potcar.psctr_parameters[-1].strip() == "local part"
    assert potcar.psctr_parameters[0].strip() == "PAW_PBE W 08Apr2002"
    assert np.isclose(potcar.k_max, 68.0301355888413)
    assert np.isclose(potcar.fourier_coefficients[0], 0.12147212e3)
    assert np.isclose(potcar.fourier_coefficients[-1], 0.23490619e0)
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

    assert potcar.psctr_parameters[-1].strip() == "local part"
    assert potcar.psctr_parameters[0].strip() == "PAW_PBE Xe 07Sep2000"
    assert np.isclose(potcar.k_max, 75.5890395431569)
    assert np.isclose(potcar.fourier_coefficients[0], 0.29945148e3)
    assert np.isclose(potcar.fourier_coefficients[-1], 0.25592291e0)
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
