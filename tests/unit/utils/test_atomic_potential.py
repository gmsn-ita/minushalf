"""
Test atomic potential
"""
import numpy as np
from minushalf.utils import Vtotal, AtomicPotential
from minushalf.softwares.vasp import Potcar


def test_occupy_potential_ag(file_path):
    """
    Test occupy potential function for Ag
    """
    vtotal_path = file_path('/Ag/VTOTAL')
    vtotal_occ_path = file_path('/Ag/VTOTAL_OCC')
    potcar_path = file_path('/Ag/POTCAR')

    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, potcar)

    cut = 2.5110905804852792
    amplitude = 1.2884975612107934
    occupy_potential = atomic_potential.occupy_potential(cut, amplitude)

    assert np.isclose(occupy_potential[0], -2.066447735049675e-05)
    assert np.isclose(occupy_potential[-1], 0.0)
    assert len(occupy_potential) == 1171


def test_correct_potential_ag(file_path):
    """
    Test correct potential function for Ag
    """
    vtotal_path = file_path('/Ag/VTOTAL')
    vtotal_occ_path = file_path('/Ag/VTOTAL_OCC')
    potcar_path = file_path('/Ag/POTCAR')

    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, potcar)

    cut = 2.5110905804852792
    amplitude = 1.2884975612107934
    is_conduction = True
    corrected_potential = atomic_potential.correct_potential(
        cut, amplitude, is_conduction)

    assert np.isclose(corrected_potential[0], 42.134427008617926)
    assert np.isclose(corrected_potential[-1], 3.3388329555573573)
    assert len(corrected_potential) == 1000


def test_occupy_potential_c(file_path):
    """
    Test occupy potential function for C
    """
    vtotal_path = file_path('/C/VTOTAL')
    vtotal_occ_path = file_path('/C/VTOTAL_OCC')
    potcar_path = file_path('/C/POTCAR')

    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, potcar)

    cut = 2.9292522566374504
    amplitude = 1.7855489612373474
    occupy_potential = atomic_potential.occupy_potential(cut, amplitude)

    assert np.isclose(occupy_potential[0], -0.0015258608114589413)
    assert np.isclose(occupy_potential[-1], 0.0)
    assert len(occupy_potential) == 1006


def test_correct_potential_c(file_path):
    """
    Test correct potential function for C
    """
    vtotal_path = file_path('/C/VTOTAL')
    vtotal_occ_path = file_path('/C/VTOTAL_OCC')
    potcar_path = file_path('/C/POTCAR')

    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, potcar)

    cut = 2.9292522566374504
    amplitude = 1.7855489612373474
    is_conduction = True
    corrected_potential = atomic_potential.correct_potential(
        cut, amplitude, is_conduction)

    assert np.isclose(corrected_potential[0], 82.71034211169965)
    assert np.isclose(corrected_potential[-1], 2.222216630996971)
    assert len(corrected_potential) == 1000


def test_occupy_potential_er(file_path):
    """
    Test occupy potential function for Er
    """
    vtotal_path = file_path('/Er/VTOTAL')
    vtotal_occ_path = file_path('/Er/VTOTAL_OCC')
    potcar_path = file_path('/Er/POTCAR')

    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, potcar)

    cut = 3.5498933827077526
    amplitude = 1.216568867569384
    occupy_potential = atomic_potential.occupy_potential(cut, amplitude)

    assert np.isclose(occupy_potential[0], -2.0867117713605294e-05)
    assert np.isclose(occupy_potential[-1], 0.0)
    assert len(occupy_potential) == 1200


def test_correct_potential_er(file_path):
    """
    Test correct potential function for Er
    """
    vtotal_path = file_path('/Er/VTOTAL')
    vtotal_occ_path = file_path('/Er/VTOTAL_OCC')
    potcar_path = file_path('/Er/POTCAR')

    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, potcar)

    cut = 3.5498933827077526
    amplitude = 1.216568867569384
    is_conduction = False
    corrected_potential = atomic_potential.correct_potential(
        cut, amplitude, is_conduction)

    assert np.isclose(corrected_potential[0], -82.04984588007534)
    assert np.isclose(corrected_potential[-1], 3.3337994605567602)
    assert len(corrected_potential) == 1000


def test_occupy_potential_f(file_path):
    """
    Test occupy potential function for F
    """
    vtotal_path = file_path('/F/VTOTAL')
    vtotal_occ_path = file_path('/F/VTOTAL_OCC')
    potcar_path = file_path('/F/POTCAR')

    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, potcar)

    cut = 3.0432090491580475
    amplitude = 1.3026752967135127
    occupy_potential = atomic_potential.occupy_potential(cut, amplitude)

    assert np.isclose(occupy_potential[0], -0.0006887447208360254)
    assert np.isclose(occupy_potential[-1], 0.0)
    assert len(occupy_potential) == 1038


def test_correct_potential_f(file_path):
    """
    Test correct potential function for F
    """
    vtotal_path = file_path('/F/VTOTAL')
    vtotal_occ_path = file_path('/F/VTOTAL_OCC')
    potcar_path = file_path('/F/POTCAR')

    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, potcar)

    cut = 3.0432090491580475
    amplitude = 1.3026752967135127
    is_conduction = True
    corrected_potential = atomic_potential.correct_potential(
        cut, amplitude, is_conduction)

    assert np.isclose(corrected_potential[0], 78.5746517531089)
    assert np.isclose(corrected_potential[-1], 1.1111107999962797)
    assert len(corrected_potential) == 1000


def test_occupy_potential_fe(file_path):
    """
    Test occupy potential function for Fe
    """
    vtotal_path = file_path('/Fe/VTOTAL')
    vtotal_occ_path = file_path('/Fe/VTOTAL_OCC')
    potcar_path = file_path('/Fe/POTCAR')

    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, potcar)

    cut = 3.8163735891237063
    amplitude = 1.3643431827480144
    occupy_potential = atomic_potential.occupy_potential(cut, amplitude)

    assert np.isclose(occupy_potential[0], -0.00010067271074148367)
    assert np.isclose(occupy_potential[-1], 0.0)
    assert len(occupy_potential) == 1123


def test_correct_potential_fe(file_path):
    """
    Test correct potential function for Fe
    """
    vtotal_path = file_path('/Fe/VTOTAL')
    vtotal_occ_path = file_path('/Fe/VTOTAL_OCC')
    potcar_path = file_path('/Fe/POTCAR')

    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, potcar)

    cut = 3.8163735891237063
    amplitude = 1.3643431827480144
    is_conduction = True
    corrected_potential = atomic_potential.correct_potential(
        cut, amplitude, is_conduction)

    assert np.isclose(corrected_potential[0], 120.74779607321234)
    assert np.isclose(corrected_potential[-1], 4.419636704868847)
    assert len(corrected_potential) == 1000


def test_occupy_potential_ga(file_path):
    """
    Test occupy potential function for Ga
    """
    vtotal_path = file_path('/Ga/VTOTAL')
    vtotal_occ_path = file_path('/Ga/VTOTAL_OCC')
    potcar_path = file_path('/Ga/POTCAR')

    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, potcar)

    cut = 2.8665035768110703
    amplitude = 1.2021549956518776
    occupy_potential = atomic_potential.occupy_potential(cut, amplitude)

    assert np.isclose(occupy_potential[0], -2.3881908615319164e-05)
    assert np.isclose(occupy_potential[-1], 0.0)
    assert len(occupy_potential) == 1137


def test_correct_potential_ga(file_path):
    """
    Test correct potential function for Ga
    """
    vtotal_path = file_path('/Ga/VTOTAL')
    vtotal_occ_path = file_path('/Ga/VTOTAL_OCC')
    potcar_path = file_path('/Ga/POTCAR')

    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, potcar)

    cut = 2.8665035768110703
    amplitude = 1.2021549956518776
    is_conduction = True
    corrected_potential = atomic_potential.correct_potential(
        cut, amplitude, is_conduction)

    assert np.isclose(corrected_potential[0], 40.22661269132331)
    assert np.isclose(corrected_potential[-1], 2.218316919876328)
    assert len(corrected_potential) == 1000


def test_occupy_potential_h(file_path):
    """
    Test occupy potential function for H
    """
    vtotal_path = file_path('/H/VTOTAL')
    vtotal_occ_path = file_path('/H/VTOTAL_OCC')
    potcar_path = file_path('/H/POTCAR')

    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, potcar)

    cut = 3.4535972492030687
    amplitude = 1.810139694953808
    occupy_potential = atomic_potential.occupy_potential(cut, amplitude)

    assert np.isclose(occupy_potential[0], 0.018199220885673502)
    assert np.isclose(occupy_potential[-1], 0.0)
    assert len(occupy_potential) == 863


def test_correct_potential_h(file_path):
    """
    Test correct potential function for H
    """
    vtotal_path = file_path('/H/VTOTAL')
    vtotal_occ_path = file_path('/H/VTOTAL_OCC')
    potcar_path = file_path('/H/POTCAR')

    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, potcar)

    cut = 3.4535972492030687
    amplitude = 1.810139694953808
    is_conduction = False
    corrected_potential = atomic_potential.correct_potential(
        cut, amplitude, is_conduction)

    assert np.isclose(corrected_potential[0], -118.03778676634806)
    assert np.isclose(corrected_potential[-1], 4.4443560211104325)
    assert len(corrected_potential) == 1000


def test_occupy_potential_he(file_path):
    """
    Test occupy potential function for He
    """
    vtotal_path = file_path('/He/VTOTAL')
    vtotal_occ_path = file_path('/He/VTOTAL_OCC')
    potcar_path = file_path('/He/POTCAR')

    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, potcar)

    cut = 3.0243056287111063
    amplitude = 1.5660092595324393
    occupy_potential = atomic_potential.occupy_potential(cut, amplitude)

    assert np.isclose(occupy_potential[0], 0.04016985857986955)
    assert np.isclose(occupy_potential[-1], 0.0)
    assert len(occupy_potential) == 918


def test_correct_potential_he(file_path):
    """
    Test correct potential function for He
    """
    vtotal_path = file_path('/He/VTOTAL')
    vtotal_occ_path = file_path('/He/VTOTAL_OCC')
    potcar_path = file_path('/He/POTCAR')

    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, potcar)

    cut = 3.0243056287111063
    amplitude = 1.5660092595324393
    is_conduction = True
    corrected_potential = atomic_potential.correct_potential(
        cut, amplitude, is_conduction)

    assert np.isclose(corrected_potential[0], 98.95795397153951)
    assert np.isclose(corrected_potential[-1], 4.4337736399639365)
    assert len(corrected_potential) == 1000


def test_occupy_potential_ir(file_path):
    """
    Test occupy potential function for Ir
    """
    vtotal_path = file_path('/Ir/VTOTAL')
    vtotal_occ_path = file_path('/Ir/VTOTAL_OCC')
    potcar_path = file_path('/Ir/POTCAR')

    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, potcar)

    cut = 3.5578895515171602
    amplitude = 1.8042252321570866
    occupy_potential = atomic_potential.occupy_potential(cut, amplitude)

    assert np.isclose(occupy_potential[0], -1.1702208195045162e-05)
    assert np.isclose(occupy_potential[-1], 0.0)
    assert len(occupy_potential) == 1210


def test_correct_potential_ir(file_path):
    """
    Test correct potential function for Ir
    """
    vtotal_path = file_path('/Ir/VTOTAL')
    vtotal_occ_path = file_path('/Ir/VTOTAL_OCC')
    potcar_path = file_path('/Ir/POTCAR')

    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, potcar)

    cut = 3.5578895515171602
    amplitude = 1.8042252321570866
    is_conduction = False
    corrected_potential = atomic_potential.correct_potential(
        cut, amplitude, is_conduction)

    assert np.isclose(corrected_potential[0], -116.10754182680598)
    assert np.isclose(corrected_potential[-1], 4.412557155252129)
    assert len(corrected_potential) == 1000


def test_occupy_potential_mg(file_path):
    """
    Test occupy potential function for Mg
    """
    vtotal_path = file_path('/Mg/VTOTAL')
    vtotal_occ_path = file_path('/Mg/VTOTAL_OCC')
    potcar_path = file_path('/Mg/POTCAR')

    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, potcar)

    cut = 2.467345936958237
    amplitude = 1.2013972481132784
    occupy_potential = atomic_potential.occupy_potential(cut, amplitude)

    assert np.isclose(occupy_potential[0], 6.962814910261916e-05)
    assert np.isclose(occupy_potential[-1], 0.0)
    assert len(occupy_potential) == 1061


def test_correct_potential_mg(file_path):
    """
    Test correct potential function for Mg
    """
    vtotal_path = file_path('/Mg/VTOTAL')
    vtotal_occ_path = file_path('/Mg/VTOTAL_OCC')
    potcar_path = file_path('/Mg/POTCAR')

    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, potcar)

    cut = 2.467345936958237
    amplitude = 1.2013972481132784
    is_conduction = True
    corrected_potential = atomic_potential.correct_potential(
        cut, amplitude, is_conduction)

    assert np.isclose(corrected_potential[0], 25.954627435540964)
    assert np.isclose(corrected_potential[-1], 3.3338022720905234)
    assert len(corrected_potential) == 1000


def test_occupy_potential_mo(file_path):
    """
    Test occupy potential function for Mo
    """
    vtotal_path = file_path('/Mo/VTOTAL')
    vtotal_occ_path = file_path('/Mo/VTOTAL_OCC')
    potcar_path = file_path('/Mo/POTCAR')

    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, potcar)

    cut = 2.180064656117178
    amplitude = 1.42058347245096
    occupy_potential = atomic_potential.occupy_potential(cut, amplitude)

    assert np.isclose(occupy_potential[0], -1.9719921147957657e-05)
    assert np.isclose(occupy_potential[-1], 0.0)
    assert len(occupy_potential) == 1162


def test_correct_potential_mo(file_path):
    """
    Test correct potential function for Mo
    """
    vtotal_path = file_path('/Mo/VTOTAL')
    vtotal_occ_path = file_path('/Mo/VTOTAL_OCC')
    potcar_path = file_path('/Mo/POTCAR')

    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, potcar)

    cut = 2.180064656117178
    amplitude = 1.42058347245096
    is_conduction = True
    corrected_potential = atomic_potential.correct_potential(
        cut, amplitude, is_conduction)

    assert np.isclose(corrected_potential[0], 25.349795031355335)
    assert np.isclose(corrected_potential[-1], 1.1292868819325126)
    assert len(corrected_potential) == 1000


def test_occupy_potential_n(file_path):
    """
    Test occupy potential function for N
    """
    vtotal_path = file_path('/N/VTOTAL')
    vtotal_occ_path = file_path('/N/VTOTAL_OCC')
    potcar_path = file_path('/N/POTCAR')

    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, potcar)

    cut = 2.225580070973389
    amplitude = 1.9618849835282237
    occupy_potential = atomic_potential.occupy_potential(cut, amplitude)

    assert np.isclose(occupy_potential[0], -0.0013927531179723342)
    assert np.isclose(occupy_potential[-1], 0.0)
    assert len(occupy_potential) == 1018


def test_correct_potential_n(file_path):
    """
    Test correct potential function for N
    """
    vtotal_path = file_path('/N/VTOTAL')
    vtotal_occ_path = file_path('/N/VTOTAL_OCC')
    potcar_path = file_path('/N/POTCAR')

    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, potcar)

    cut = 2.225580070973389
    amplitude = 1.9618849835282237
    is_conduction = True
    corrected_potential = atomic_potential.correct_potential(
        cut, amplitude, is_conduction)

    assert np.isclose(corrected_potential[0], 50.754467050795675)
    assert np.isclose(corrected_potential[-1], 2.2232765357932323)
    assert len(corrected_potential) == 1000


def test_occupy_potential_ne(file_path):
    """
    Test occupy potential function for Ne
    """
    vtotal_path = file_path('/Ne/VTOTAL')
    vtotal_occ_path = file_path('/Ne/VTOTAL_OCC')
    potcar_path = file_path('/Ne/POTCAR')

    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, potcar)

    cut = 2.428368581472277
    amplitude = 1.9549679403647553
    occupy_potential = atomic_potential.occupy_potential(cut, amplitude)

    assert np.isclose(occupy_potential[0], -0.0009160652192273718)
    assert np.isclose(occupy_potential[-1], 0.0)
    assert len(occupy_potential) == 1047


def test_correct_potential_ne(file_path):
    """
    Test correct potential function for Ne
    """
    vtotal_path = file_path('/Ne/VTOTAL')
    vtotal_occ_path = file_path('/Ne/VTOTAL_OCC')
    potcar_path = file_path('/Ne/POTCAR')

    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, potcar)

    cut = 2.428368581472277
    amplitude = 1.9549679403647553
    is_conduction = True
    corrected_potential = atomic_potential.correct_potential(
        cut, amplitude, is_conduction)

    assert np.isclose(corrected_potential[0], 71.53985054646311)
    assert np.isclose(corrected_potential[-1], 4.444443311513178)
    assert len(corrected_potential) == 1000


def test_occupy_potential_ru(file_path):
    """
    Test occupy potential function for Ru
    """
    vtotal_path = file_path('/Ru/VTOTAL')
    vtotal_occ_path = file_path('/Ru/VTOTAL_OCC')
    potcar_path = file_path('/Ru/POTCAR')

    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, potcar)

    cut = 3.7008862737009958
    amplitude = 1.6105876308829457
    occupy_potential = atomic_potential.occupy_potential(cut, amplitude)

    assert np.isclose(occupy_potential[0], -2.4087640516951823e-05)
    assert np.isclose(occupy_potential[-1], 0.0)
    assert len(occupy_potential) == 1165


def test_correct_potential_ru(file_path):
    """
    Test correct potential function for Ru
    """
    vtotal_path = file_path('/Ru/VTOTAL')
    vtotal_occ_path = file_path('/Ru/VTOTAL_OCC')
    potcar_path = file_path('/Ru/POTCAR')

    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, potcar)

    cut = 3.7008862737009958
    amplitude = 1.6105876308829457
    is_conduction = False
    corrected_potential = atomic_potential.correct_potential(
        cut, amplitude, is_conduction)

    assert np.isclose(corrected_potential[0], -114.60821673075468)
    assert np.isclose(corrected_potential[-1], 4.499762390411456)
    assert len(corrected_potential) == 1000


def test_occupy_potential_s(file_path):
    """
    Test occupy potential function for S
    """
    vtotal_path = file_path('/S/VTOTAL')
    vtotal_occ_path = file_path('/S/VTOTAL_OCC')
    potcar_path = file_path('/S/POTCAR')

    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, potcar)

    cut = 2.704287004646292
    amplitude = 1.5524108053068066
    occupy_potential = atomic_potential.occupy_potential(cut, amplitude)

    assert np.isclose(occupy_potential[0], -0.0001294566715730413)
    assert np.isclose(occupy_potential[-1], 0.0)
    assert len(occupy_potential) == 1084


def test_correct_potential_s(file_path):
    """
    Test correct potential function for S
    """
    vtotal_path = file_path('/S/VTOTAL')
    vtotal_occ_path = file_path('/S/VTOTAL_OCC')
    potcar_path = file_path('/S/POTCAR')

    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, potcar)

    cut = 2.704287004646292
    amplitude = 1.5524108053068066
    is_conduction = False
    corrected_potential = atomic_potential.correct_potential(
        cut, amplitude, is_conduction)

    assert np.isclose(corrected_potential[0], -47.430987087457154)
    assert np.isclose(corrected_potential[-1], 4.445449320860461)
    assert len(corrected_potential) == 1000


def test_occupy_potential_sn(file_path):
    """
    Test occupy potential function for Sn
    """
    vtotal_path = file_path('/Sn/VTOTAL')
    vtotal_occ_path = file_path('/Sn/VTOTAL_OCC')
    potcar_path = file_path('/Sn/POTCAR')

    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, potcar)

    cut = 3.5439491123023767
    amplitude = 1.5318524949299626
    occupy_potential = atomic_potential.occupy_potential(cut, amplitude)

    assert np.isclose(occupy_potential[0], -1.2935683695324023e-05)
    assert np.isclose(occupy_potential[-1], 0.0)
    assert len(occupy_potential) == 1175


def test_correct_potential_sn(file_path):
    """
    Test correct potential function for Sn
    """
    vtotal_path = file_path('/Sn/VTOTAL')
    vtotal_occ_path = file_path('/Sn/VTOTAL_OCC')
    potcar_path = file_path('/Sn/POTCAR')

    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, potcar)

    cut = 3.5439491123023767
    amplitude = 1.5318524949299626
    is_conduction = False
    corrected_potential = atomic_potential.correct_potential(
        cut, amplitude, is_conduction)

    assert np.isclose(corrected_potential[0], -80.53606842010086)
    assert np.isclose(corrected_potential[-1], 1.4141791515446953)
    assert len(corrected_potential) == 1000


def test_occupy_potential_tc(file_path):
    """
    Test occupy potential function for Tc
    """
    vtotal_path = file_path('/Tc/VTOTAL')
    vtotal_occ_path = file_path('/Tc/VTOTAL_OCC')
    potcar_path = file_path('/Tc/POTCAR')

    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, potcar)

    cut = 3.761467176330928
    amplitude = 1.2293238434175695
    occupy_potential = atomic_potential.occupy_potential(cut, amplitude)

    assert np.isclose(occupy_potential[0], -1.7781298095659134e-05)
    assert np.isclose(occupy_potential[-1], 0.0)
    assert len(occupy_potential) == 1163


def test_correct_potential_tc(file_path):
    """
    Test correct potential function for Tc
    """
    vtotal_path = file_path('/Tc/VTOTAL')
    vtotal_occ_path = file_path('/Tc/VTOTAL_OCC')
    potcar_path = file_path('/Tc/POTCAR')

    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, potcar)

    cut = 3.761467176330928
    amplitude = 1.2293238434175695
    is_conduction = True
    corrected_potential = atomic_potential.correct_potential(
        cut, amplitude, is_conduction)

    assert np.isclose(corrected_potential[0], 93.21778221959654)
    assert np.isclose(corrected_potential[-1], 1.2892575561440929)
    assert len(corrected_potential) == 1000


def test_occupy_potential_v(file_path):
    """
    Test occupy potential function for V
    """
    vtotal_path = file_path('/V/VTOTAL')
    vtotal_occ_path = file_path('/V/VTOTAL_OCC')
    potcar_path = file_path('/V/POTCAR')

    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, potcar)

    cut = 3.0666023646202416
    amplitude = 1.9240605455394564
    occupy_potential = atomic_potential.occupy_potential(cut, amplitude)

    assert np.isclose(occupy_potential[0], -0.00015413596312796796)
    assert np.isclose(occupy_potential[-1], 0.0)
    assert len(occupy_potential) == 1113


def test_correct_potential_v(file_path):
    """
    Test correct potential function for V
    """
    vtotal_path = file_path('/V/VTOTAL')
    vtotal_occ_path = file_path('/V/VTOTAL_OCC')
    potcar_path = file_path('/V/POTCAR')

    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, potcar)

    cut = 3.0666023646202416
    amplitude = 1.9240605455394564
    is_conduction = True
    corrected_potential = atomic_potential.correct_potential(
        cut, amplitude, is_conduction)

    assert np.isclose(corrected_potential[0], 92.55058321386603)
    assert np.isclose(corrected_potential[-1], 2.221542666741145)
    assert len(corrected_potential) == 1000


def test_occupy_potential_w(file_path):
    """
    Test occupy potential function for W
    """
    vtotal_path = file_path('/W/VTOTAL')
    vtotal_occ_path = file_path('/W/VTOTAL_OCC')
    potcar_path = file_path('/W/POTCAR')

    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, potcar)

    cut = 2.4557084902787274
    amplitude = 1.0063540842881846
    occupy_potential = atomic_potential.occupy_potential(cut, amplitude)

    assert np.isclose(occupy_potential[0], -8.210012297487316e-06)
    assert np.isclose(occupy_potential[-1], 0.0)
    assert len(occupy_potential) == 1207


def test_correct_potential_w(file_path):
    """
    Test correct potential function for W
    """
    vtotal_path = file_path('/W/VTOTAL')
    vtotal_occ_path = file_path('/W/VTOTAL_OCC')
    potcar_path = file_path('/W/POTCAR')

    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, potcar)

    cut = 2.4557084902787274
    amplitude = 1.0063540842881846
    is_conduction = True
    corrected_potential = atomic_potential.correct_potential(
        cut, amplitude, is_conduction)

    assert np.isclose(corrected_potential[0], 28.22421669671955)
    assert np.isclose(corrected_potential[-1], 4.444104106630428)
    assert len(corrected_potential) == 1000


def test_occupy_potential_xe(file_path):
    """
    Test occupy potential function for Xe
    """
    vtotal_path = file_path('/Xe/VTOTAL')
    vtotal_occ_path = file_path('/Xe/VTOTAL_OCC')
    potcar_path = file_path('/Xe/POTCAR')

    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, potcar)

    cut = 2.9783715913046125
    amplitude = 1.0428254230615635
    occupy_potential = atomic_potential.occupy_potential(cut, amplitude)

    assert np.isclose(occupy_potential[0], -9.194496883803786e-06)
    assert np.isclose(occupy_potential[-1], 0.0)
    assert len(occupy_potential) == 1182


def test_correct_potential_xe(file_path):
    """
    Test correct potential function for Xe
    """
    vtotal_path = file_path('/Xe/VTOTAL')
    vtotal_occ_path = file_path('/Xe/VTOTAL_OCC')
    potcar_path = file_path('/Xe/POTCAR')

    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, potcar)

    cut = 2.9783715913046125
    amplitude = 1.0428254230615635
    is_conduction = True
    corrected_potential = atomic_potential.correct_potential(
        cut, amplitude, is_conduction)

    assert np.isclose(corrected_potential[0], 43.33294744926906)
    assert np.isclose(corrected_potential[-1], 6.13351807697526)
    assert len(corrected_potential) == 1000


def test_get_lines_corrected_potential_xe(file_path):
    """
    Test get lines of corrected potential function for Xe
    """
    vtotal_path = file_path('/Xe/VTOTAL')
    vtotal_occ_path = file_path('/Xe/VTOTAL_OCC')
    potcar_path = file_path('/Xe/POTCAR')

    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, potcar)

    cut = 2.9783715913046125
    amplitude = 1.0428254230615635
    is_conduction = True
    corrected_potential = atomic_potential.correct_potential(
        cut, amplitude, is_conduction)
    lines = atomic_potential.get_corrected_file_lines(corrected_potential)

    assert lines[3].strip() == "PAW_PBE Xe 21Sep2222"
    assert lines[-1].strip() == "End of Dataset"
    assert len(lines) == 2040


def test_get_lines_corrected_potential_w(file_path):
    """
    Test get lines of corrected potential function for W
    """
    vtotal_path = file_path('/W/VTOTAL')
    vtotal_occ_path = file_path('/W/VTOTAL_OCC')
    potcar_path = file_path('/W/POTCAR')

    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, potcar)

    cut = 2.4557084902787274
    amplitude = 1.0063540842881846
    is_conduction = True
    corrected_potential = atomic_potential.correct_potential(
        cut, amplitude, is_conduction)
    lines = atomic_potential.get_corrected_file_lines(corrected_potential)

    assert lines[3].strip() == "PAW_PBE W 41Apr4444"
    assert lines[-1].strip() == "End of Dataset"
    assert len(lines) == 2586
