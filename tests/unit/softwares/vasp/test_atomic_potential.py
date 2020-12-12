"""
Test atomic potential
"""
import numpy as np
from minushalf.atomic import InputFile, Vtotal
from minushalf.softwares.vasp import AtomicPotential, Potcar


def test_occupy_potential_ag(file_path):
    """
    Test occupy potential function for Ag
    """
    inp_path = file_path('/Ag/INP_OCC')
    vtotal_path = file_path('/Ag/VTOTAL')
    vtotal_occ_path = file_path('/Ag/VTOTAL_OCC')
    potcar_path = file_path('/Ag/POTCAR')

    input_file = InputFile.from_file(inp_path)
    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, input_file, potcar)

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
    inp_path = file_path('/Ag/INP_OCC')
    vtotal_path = file_path('/Ag/VTOTAL')
    vtotal_occ_path = file_path('/Ag/VTOTAL_OCC')
    potcar_path = file_path('/Ag/POTCAR')

    input_file = InputFile.from_file(inp_path)
    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, input_file, potcar)

    cut = 2.5110905804852792
    amplitude = 1.2884975612107934
    is_conduction = True
    corrected_potential = atomic_potential.correct_potential(
        cut, amplitude, is_conduction)

    assert np.isclose(corrected_potential[0], 115.01526160861792)
    assert np.isclose(corrected_potential[-1], 0.3483004586739098)
    assert len(corrected_potential) == 1000


def test_occupy_potential_c(file_path):
    """
    Test occupy potential function for C
    """
    inp_path = file_path('/C/INP_OCC')
    vtotal_path = file_path('/C/VTOTAL')
    vtotal_occ_path = file_path('/C/VTOTAL_OCC')
    potcar_path = file_path('/C/POTCAR')

    input_file = InputFile.from_file(inp_path)
    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, input_file, potcar)

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
    inp_path = file_path('/C/INP_OCC')
    vtotal_path = file_path('/C/VTOTAL')
    vtotal_occ_path = file_path('/C/VTOTAL_OCC')
    potcar_path = file_path('/C/POTCAR')

    input_file = InputFile.from_file(inp_path)
    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, input_file, potcar)

    cut = 2.9292522566374504
    amplitude = 1.7855489612373474
    is_conduction = True
    corrected_potential = atomic_potential.correct_potential(
        cut, amplitude, is_conduction)

    assert np.isclose(corrected_potential[0], 98.80257671169966)
    assert np.isclose(corrected_potential[-1], 0.0469647743921682)
    assert len(corrected_potential) == 1000


def test_occupy_potential_er(file_path):
    """
    Test occupy potential function for Er
    """
    inp_path = file_path('/Er/INP_OCC')
    vtotal_path = file_path('/Er/VTOTAL')
    vtotal_occ_path = file_path('/Er/VTOTAL_OCC')
    potcar_path = file_path('/Er/POTCAR')

    input_file = InputFile.from_file(inp_path)
    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, input_file, potcar)

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
    inp_path = file_path('/Er/INP_OCC')
    vtotal_path = file_path('/Er/VTOTAL')
    vtotal_occ_path = file_path('/Er/VTOTAL_OCC')
    potcar_path = file_path('/Er/POTCAR')

    input_file = InputFile.from_file(inp_path)
    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, input_file, potcar)

    cut = 3.5498933827077526
    amplitude = 1.216568867569384
    is_conduction = False
    corrected_potential = atomic_potential.correct_potential(
        cut, amplitude, is_conduction)

    assert np.isclose(corrected_potential[0], 247.34325191992465)
    assert np.isclose(corrected_potential[-1], 0.7748708762400962)
    assert len(corrected_potential) == 1000


def test_occupy_potential_f(file_path):
    """
    Test occupy potential function for F
    """
    inp_path = file_path('/F/INP_OCC')
    vtotal_path = file_path('/F/VTOTAL')
    vtotal_occ_path = file_path('/F/VTOTAL_OCC')
    potcar_path = file_path('/F/POTCAR')

    input_file = InputFile.from_file(inp_path)
    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, input_file, potcar)

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
    inp_path = file_path('/F/INP_OCC')
    vtotal_path = file_path('/F/VTOTAL')
    vtotal_occ_path = file_path('/F/VTOTAL_OCC')
    potcar_path = file_path('/F/POTCAR')

    input_file = InputFile.from_file(inp_path)
    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, input_file, potcar)

    cut = 3.0432090491580475
    amplitude = 1.3026752967135127
    is_conduction = True
    corrected_potential = atomic_potential.correct_potential(
        cut, amplitude, is_conduction)

    assert np.isclose(corrected_potential[0], 100.3215204531089)
    assert np.isclose(corrected_potential[-1], 0.08422915502323297)
    assert len(corrected_potential) == 1000


def test_occupy_potential_fe(file_path):
    """
    Test occupy potential function for Fe
    """
    inp_path = file_path('/Fe/INP_OCC')
    vtotal_path = file_path('/Fe/VTOTAL')
    vtotal_occ_path = file_path('/Fe/VTOTAL_OCC')
    potcar_path = file_path('/Fe/POTCAR')

    input_file = InputFile.from_file(inp_path)
    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, input_file, potcar)

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
    inp_path = file_path('/Fe/INP_OCC')
    vtotal_path = file_path('/Fe/VTOTAL')
    vtotal_occ_path = file_path('/Fe/VTOTAL_OCC')
    potcar_path = file_path('/Fe/POTCAR')

    input_file = InputFile.from_file(inp_path)
    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, input_file, potcar)

    cut = 3.8163735891237063
    amplitude = 1.3643431827480144
    is_conduction = True
    corrected_potential = atomic_potential.correct_potential(
        cut, amplitude, is_conduction)

    assert np.isclose(corrected_potential[0], 202.79657167321233)
    assert np.isclose(corrected_potential[-1], 0.2200274503723977)
    assert len(corrected_potential) == 1000


def test_occupy_potential_ga(file_path):
    """
    Test occupy potential function for Ga
    """
    inp_path = file_path('/Ga/INP_OCC')
    vtotal_path = file_path('/Ga/VTOTAL')
    vtotal_occ_path = file_path('/Ga/VTOTAL_OCC')
    potcar_path = file_path('/Ga/POTCAR')

    input_file = InputFile.from_file(inp_path)
    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, input_file, potcar)

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
    inp_path = file_path('/Ga/INP_OCC')
    vtotal_path = file_path('/Ga/VTOTAL')
    vtotal_occ_path = file_path('/Ga/VTOTAL_OCC')
    potcar_path = file_path('/Ga/POTCAR')

    input_file = InputFile.from_file(inp_path)
    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, input_file, potcar)

    cut = 2.8665035768110703
    amplitude = 1.2021549956518776
    is_conduction = True
    corrected_potential = atomic_potential.correct_potential(
        cut, amplitude, is_conduction)

    assert np.isclose(corrected_potential[0], 95.2423412913233)
    assert np.isclose(corrected_potential[-1], 0.10633063376322688)
    assert len(corrected_potential) == 1000


def test_occupy_potential_h(file_path):
    """
    Test occupy potential function for H
    """
    inp_path = file_path('/H/INP_OCC')
    vtotal_path = file_path('/H/VTOTAL')
    vtotal_occ_path = file_path('/H/VTOTAL_OCC')
    potcar_path = file_path('/H/POTCAR')

    input_file = InputFile.from_file(inp_path)
    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, input_file, potcar)

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
    inp_path = file_path('/H/INP_OCC')
    vtotal_path = file_path('/H/VTOTAL')
    vtotal_occ_path = file_path('/H/VTOTAL_OCC')
    potcar_path = file_path('/H/POTCAR')

    input_file = InputFile.from_file(inp_path)
    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, input_file, potcar)

    cut = 3.4535972492030687
    amplitude = 1.810139694953808
    is_conduction = False
    corrected_potential = atomic_potential.correct_potential(
        cut, amplitude, is_conduction)

    assert np.isclose(corrected_potential[0], -118.93495256634807)
    assert np.isclose(corrected_potential[-1], 0.00622937466765834)
    assert len(corrected_potential) == 1000


def test_occupy_potential_he(file_path):
    """
    Test occupy potential function for He
    """
    inp_path = file_path('/He/INP_OCC')
    vtotal_path = file_path('/He/VTOTAL')
    vtotal_occ_path = file_path('/He/VTOTAL_OCC')
    potcar_path = file_path('/He/POTCAR')

    input_file = InputFile.from_file(inp_path)
    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, input_file, potcar)

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
    inp_path = file_path('/He/INP_OCC')
    vtotal_path = file_path('/He/VTOTAL')
    vtotal_occ_path = file_path('/He/VTOTAL_OCC')
    potcar_path = file_path('/He/POTCAR')

    input_file = InputFile.from_file(inp_path)
    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, input_file, potcar)

    cut = 3.0243056287111063
    amplitude = 1.5660092595324393
    is_conduction = True
    corrected_potential = atomic_potential.correct_potential(
        cut, amplitude, is_conduction)

    assert np.isclose(corrected_potential[0], 100.4524327715395)
    assert np.isclose(corrected_potential[-1], 0.012553526993934958)
    assert len(corrected_potential) == 1000


def test_occupy_potential_ir(file_path):
    """
    Test occupy potential function for Ir
    """
    inp_path = file_path('/Ir/INP_OCC')
    vtotal_path = file_path('/Ir/VTOTAL')
    vtotal_occ_path = file_path('/Ir/VTOTAL_OCC')
    potcar_path = file_path('/Ir/POTCAR')

    input_file = InputFile.from_file(inp_path)
    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, input_file, potcar)

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
    inp_path = file_path('/Ir/INP_OCC')
    vtotal_path = file_path('/Ir/VTOTAL')
    vtotal_occ_path = file_path('/Ir/VTOTAL_OCC')
    potcar_path = file_path('/Ir/POTCAR')

    input_file = InputFile.from_file(inp_path)
    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, input_file, potcar)

    cut = 3.5578895515171602
    amplitude = 1.8042252321570866
    is_conduction = False
    corrected_potential = atomic_potential.correct_potential(
        cut, amplitude, is_conduction)

    assert np.isclose(corrected_potential[0], 6.140173773194022)
    assert np.isclose(corrected_potential[-1], 0.3162435044942796)
    assert len(corrected_potential) == 1000


def test_occupy_potential_mg(file_path):
    """
    Test occupy potential function for Mg
    """
    inp_path = file_path('/Mg/INP_OCC')
    vtotal_path = file_path('/Mg/VTOTAL')
    vtotal_occ_path = file_path('/Mg/VTOTAL_OCC')
    potcar_path = file_path('/Mg/POTCAR')

    input_file = InputFile.from_file(inp_path)
    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, input_file, potcar)

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
    inp_path = file_path('/Mg/INP_OCC')
    vtotal_path = file_path('/Mg/VTOTAL')
    vtotal_occ_path = file_path('/Mg/VTOTAL_OCC')
    potcar_path = file_path('/Mg/POTCAR')

    input_file = InputFile.from_file(inp_path)
    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, input_file, potcar)

    cut = 2.467345936958237
    amplitude = 1.2013972481132784
    is_conduction = True
    corrected_potential = atomic_potential.correct_potential(
        cut, amplitude, is_conduction)

    assert np.isclose(corrected_potential[0], 47.34636933554096)
    assert np.isclose(corrected_potential[-1], 0.04059877663883569)
    assert len(corrected_potential) == 1000


def test_occupy_potential_mo(file_path):
    """
    Test occupy potential function for Mo
    """
    inp_path = file_path('/Mo/INP_OCC')
    vtotal_path = file_path('/Mo/VTOTAL')
    vtotal_occ_path = file_path('/Mo/VTOTAL_OCC')
    potcar_path = file_path('/Mo/POTCAR')

    input_file = InputFile.from_file(inp_path)
    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, input_file, potcar)

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
    inp_path = file_path('/Mo/INP_OCC')
    vtotal_path = file_path('/Mo/VTOTAL')
    vtotal_occ_path = file_path('/Mo/VTOTAL_OCC')
    potcar_path = file_path('/Mo/POTCAR')

    input_file = InputFile.from_file(inp_path)
    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, input_file, potcar)

    cut = 2.180064656117178
    amplitude = 1.42058347245096
    is_conduction = True
    corrected_potential = atomic_potential.correct_potential(
        cut, amplitude, is_conduction)

    assert np.isclose(corrected_potential[0], 137.53062283135532)
    assert np.isclose(corrected_potential[-1], 0.2352508564448766)
    assert len(corrected_potential) == 1000


def test_occupy_potential_n(file_path):
    """
    Test occupy potential function for N
    """
    inp_path = file_path('/N/INP_OCC')
    vtotal_path = file_path('/N/VTOTAL')
    vtotal_occ_path = file_path('/N/VTOTAL_OCC')
    potcar_path = file_path('/N/POTCAR')

    input_file = InputFile.from_file(inp_path)
    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, input_file, potcar)

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
    inp_path = file_path('/N/INP_OCC')
    vtotal_path = file_path('/N/VTOTAL')
    vtotal_occ_path = file_path('/N/VTOTAL_OCC')
    potcar_path = file_path('/N/POTCAR')

    input_file = InputFile.from_file(inp_path)
    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, input_file, potcar)

    cut = 2.225580070973389
    amplitude = 1.9618849835282237
    is_conduction = True
    corrected_potential = atomic_potential.correct_potential(
        cut, amplitude, is_conduction)

    assert np.isclose(corrected_potential[0], 67.96921965079568)
    assert np.isclose(corrected_potential[-1], 0.05861094487981112)
    assert len(corrected_potential) == 1000


def test_occupy_potential_ne(file_path):
    """
    Test occupy potential function for Ne
    """
    inp_path = file_path('/Ne/INP_OCC')
    vtotal_path = file_path('/Ne/VTOTAL')
    vtotal_occ_path = file_path('/Ne/VTOTAL_OCC')
    potcar_path = file_path('/Ne/POTCAR')

    input_file = InputFile.from_file(inp_path)
    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, input_file, potcar)

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
    inp_path = file_path('/Ne/INP_OCC')
    vtotal_path = file_path('/Ne/VTOTAL')
    vtotal_occ_path = file_path('/Ne/VTOTAL_OCC')
    potcar_path = file_path('/Ne/POTCAR')

    input_file = InputFile.from_file(inp_path)
    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, input_file, potcar)

    cut = 2.428368581472277
    amplitude = 1.9549679403647553
    is_conduction = True
    corrected_potential = atomic_potential.correct_potential(
        cut, amplitude, is_conduction)

    assert np.isclose(corrected_potential[0], 87.64705324646312)
    assert np.isclose(corrected_potential[-1], 0.12060586248801818)
    assert len(corrected_potential) == 1000


def test_occupy_potential_ru(file_path):
    """
    Test occupy potential function for Ru
    """
    inp_path = file_path('/Ru/INP_OCC')
    vtotal_path = file_path('/Ru/VTOTAL')
    vtotal_occ_path = file_path('/Ru/VTOTAL_OCC')
    potcar_path = file_path('/Ru/POTCAR')

    input_file = InputFile.from_file(inp_path)
    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, input_file, potcar)

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
    inp_path = file_path('/Ru/INP_OCC')
    vtotal_path = file_path('/Ru/VTOTAL')
    vtotal_occ_path = file_path('/Ru/VTOTAL_OCC')
    potcar_path = file_path('/Ru/POTCAR')

    input_file = InputFile.from_file(inp_path)
    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, input_file, potcar)

    cut = 3.7008862737009958
    amplitude = 1.6105876308829457
    is_conduction = False
    corrected_potential = atomic_potential.correct_potential(
        cut, amplitude, is_conduction)

    assert np.isclose(corrected_potential[0], -60.589230030754685)
    assert np.isclose(corrected_potential[-1], 0.295579883132542)
    assert len(corrected_potential) == 1000


def test_occupy_potential_s(file_path):
    """
    Test occupy potential function for S
    """
    inp_path = file_path('/S/INP_OCC')
    vtotal_path = file_path('/S/VTOTAL')
    vtotal_occ_path = file_path('/S/VTOTAL_OCC')
    potcar_path = file_path('/S/POTCAR')

    input_file = InputFile.from_file(inp_path)
    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, input_file, potcar)

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
    inp_path = file_path('/S/INP_OCC')
    vtotal_path = file_path('/S/VTOTAL')
    vtotal_occ_path = file_path('/S/VTOTAL_OCC')
    potcar_path = file_path('/S/POTCAR')

    input_file = InputFile.from_file(inp_path)
    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, input_file, potcar)

    cut = 2.704287004646292
    amplitude = 1.5524108053068066
    is_conduction = False
    corrected_potential = atomic_potential.correct_potential(
        cut, amplitude, is_conduction)

    assert np.isclose(corrected_potential[0], -26.717791487457156)
    assert np.isclose(corrected_potential[-1], 0.113274749072086)
    assert len(corrected_potential) == 1000


def test_occupy_potential_sn(file_path):
    """
    Test occupy potential function for Sn
    """
    inp_path = file_path('/Sn/INP_OCC')
    vtotal_path = file_path('/Sn/VTOTAL')
    vtotal_occ_path = file_path('/Sn/VTOTAL_OCC')
    potcar_path = file_path('/Sn/POTCAR')

    input_file = InputFile.from_file(inp_path)
    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, input_file, potcar)

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
    inp_path = file_path('/Sn/INP_OCC')
    vtotal_path = file_path('/Sn/VTOTAL')
    vtotal_occ_path = file_path('/Sn/VTOTAL_OCC')
    potcar_path = file_path('/Sn/POTCAR')

    input_file = InputFile.from_file(inp_path)
    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, input_file, potcar)

    cut = 3.5439491123023767
    amplitude = 1.5318524949299626
    is_conduction = False
    corrected_potential = atomic_potential.correct_potential(
        cut, amplitude, is_conduction)

    assert np.isclose(corrected_potential[0], -50.37416152010086)
    assert np.isclose(corrected_potential[-1], 0.18732773744957626)
    assert len(corrected_potential) == 1000


def test_occupy_potential_tc(file_path):
    """
    Test occupy potential function for Tc
    """
    inp_path = file_path('/Tc/INP_OCC')
    vtotal_path = file_path('/Tc/VTOTAL')
    vtotal_occ_path = file_path('/Tc/VTOTAL_OCC')
    potcar_path = file_path('/Tc/POTCAR')

    input_file = InputFile.from_file(inp_path)
    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, input_file, potcar)

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
    inp_path = file_path('/Tc/INP_OCC')
    vtotal_path = file_path('/Tc/VTOTAL')
    vtotal_occ_path = file_path('/Tc/VTOTAL_OCC')
    potcar_path = file_path('/Tc/POTCAR')

    input_file = InputFile.from_file(inp_path)
    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, input_file, potcar)

    cut = 3.761467176330928
    amplitude = 1.2293238434175695
    is_conduction = True
    corrected_potential = atomic_potential.correct_potential(
        cut, amplitude, is_conduction)

    assert np.isclose(corrected_potential[0], 122.29978601959654)
    assert np.isclose(corrected_potential[-1], 0.29092379877065466)
    assert len(corrected_potential) == 1000


def test_occupy_potential_v(file_path):
    """
    Test occupy potential function for V
    """
    inp_path = file_path('/V/INP_OCC')
    vtotal_path = file_path('/V/VTOTAL')
    vtotal_occ_path = file_path('/V/VTOTAL_OCC')
    potcar_path = file_path('/V/POTCAR')

    input_file = InputFile.from_file(inp_path)
    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, input_file, potcar)

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
    inp_path = file_path('/V/INP_OCC')
    vtotal_path = file_path('/V/VTOTAL')
    vtotal_occ_path = file_path('/V/VTOTAL_OCC')
    potcar_path = file_path('/V/POTCAR')

    input_file = InputFile.from_file(inp_path)
    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, input_file, potcar)

    cut = 3.0666023646202416
    amplitude = 1.9240605455394564
    is_conduction = True
    corrected_potential = atomic_potential.correct_potential(
        cut, amplitude, is_conduction)

    assert np.isclose(corrected_potential[0], 153.43341811386603)
    assert np.isclose(corrected_potential[-1], 0.18561953575232004)
    assert len(corrected_potential) == 1000


def test_occupy_potential_w(file_path):
    """
    Test occupy potential function for W
    """
    inp_path = file_path('/W/INP_OCC')
    vtotal_path = file_path('/W/VTOTAL')
    vtotal_occ_path = file_path('/W/VTOTAL_OCC')
    potcar_path = file_path('/W/POTCAR')

    input_file = InputFile.from_file(inp_path)
    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, input_file, potcar)

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
    inp_path = file_path('/W/INP_OCC')
    vtotal_path = file_path('/W/VTOTAL')
    vtotal_occ_path = file_path('/W/VTOTAL_OCC')
    potcar_path = file_path('/W/POTCAR')

    input_file = InputFile.from_file(inp_path)
    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, input_file, potcar)

    cut = 2.4557084902787274
    amplitude = 1.0063540842881846
    is_conduction = True
    corrected_potential = atomic_potential.correct_potential(
        cut, amplitude, is_conduction)

    assert np.isclose(corrected_potential[0], 145.25189229671955)
    assert np.isclose(corrected_potential[-1], 0.23495231083197574)
    assert len(corrected_potential) == 1000


def test_occupy_potential_xe(file_path):
    """
    Test occupy potential function for Xe
    """
    inp_path = file_path('/Xe/INP_OCC')
    vtotal_path = file_path('/Xe/VTOTAL')
    vtotal_occ_path = file_path('/Xe/VTOTAL_OCC')
    potcar_path = file_path('/Xe/POTCAR')

    input_file = InputFile.from_file(inp_path)
    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, input_file, potcar)

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
    inp_path = file_path('/Xe/INP_OCC')
    vtotal_path = file_path('/Xe/VTOTAL')
    vtotal_occ_path = file_path('/Xe/VTOTAL_OCC')
    potcar_path = file_path('/Xe/POTCAR')

    input_file = InputFile.from_file(inp_path)
    vtotal = Vtotal.from_file(vtotal_path)
    vtotal_occ = Vtotal.from_file(vtotal_occ_path)
    potcar = Potcar(potcar_path)
    atomic_potential = AtomicPotential(vtotal, vtotal_occ, input_file, potcar)

    cut = 2.9783715913046125
    amplitude = 1.0428254230615635
    is_conduction = True
    corrected_potential = atomic_potential.correct_potential(
        cut, amplitude, is_conduction)

    assert np.isclose(corrected_potential[0], 341.67331634926904)
    assert np.isclose(corrected_potential[-1], 0.25592177527226323)
    assert len(corrected_potential) == 1000
