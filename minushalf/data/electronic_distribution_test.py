"""
Test eletronic distribution database
"""
from collections import Counter
from minushalf.data.electronic_distribution import ElectronicDistribution


def test_electronic_distribution_h(file_path):
    """
    Test eletronic distribution of the element H
    """
    distribution_path = file_path("/H/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["H"].value)


def test_electronic_distribution_he(file_path):
    """
    Test eletronic distribution of the element He
    """
    distribution_path = file_path("/He/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["He"].value)


def test_electronic_distribution_li(file_path):
    """
    Test eletronic distribution of the element Li
    """
    distribution_path = file_path("/Li/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Li"].value)


def test_electronic_distribution_la(file_path):
    """
    Test eletronic distribution of the element La
    """
    distribution_path = file_path("/La/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["La"].value)


def test_electronic_distribution_pr(file_path):
    """
    Test eletronic distribution of the element Pr
    """
    distribution_path = file_path("/Pr/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Pr"].value)


def test_electronic_distribution_eu(file_path):
    """
    Test eletronic distribution of the element Eu
    """
    distribution_path = file_path("/Eu/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Eu"].value)


def test_electronic_distribution_tb(file_path):
    """
    Test eletronic distribution of the element Tb
    """
    distribution_path = file_path("/Tb/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Tb"].value)


def test_electronic_distribution_dy(file_path):
    """
    Test eletronic distribution of the element Dy
    """
    distribution_path = file_path("/Dy/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Dy"].value)


def test_electronic_distribution_ho(file_path):
    """
    Test eletronic distribution of the element Ho
    """
    distribution_path = file_path("Ho/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Ho"].value)


def test_electronic_distribution_fr(file_path):
    """
    Test eletronic distribution of the element Fr
    """
    distribution_path = file_path("Fr/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Fr"].value)


def test_electronic_distribution_be(file_path):
    """
    Test eletronic distribution of the element Be
    """
    distribution_path = file_path("/Be/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Be"].value)


def test_electronic_distribution_b(file_path):
    """
    Test eletronic distribution of the element B
    """
    distribution_path = file_path("/B/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["B"].value)


def test_electronic_distribution_c(file_path):
    """
    Test eletronic distribution of the element C
    """
    distribution_path = file_path("/C/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["C"].value)


def test_electronic_distribution_n(file_path):
    """
    Test eletronic distribution of the element N
    """
    distribution_path = file_path("/N/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["N"].value)


def test_electronic_distribution_o(file_path):
    """
    Test eletronic distribution of the element O
    """
    distribution_path = file_path("/O/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["O"].value)


def test_electronic_distribution_f(file_path):
    """
    Test eletronic distribution of the element F
    """
    distribution_path = file_path("/F/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["F"].value)


def test_electronic_distribution_ne(file_path):
    """
    Test eletronic distribution of the element Ne
    """
    distribution_path = file_path("/Ne/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Ne"].value)


def test_electronic_distribution_na(file_path):
    """
    Test eletronic distribution of the element Na
    """
    distribution_path = file_path("/Na/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Na"].value)


def test_electronic_distribution_mg(file_path):
    """
    Test eletronic distribution of the element Mg
    """
    distribution_path = file_path("/Mg/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Mg"].value)


def test_electronic_distribution_al(file_path):
    """
    Test eletronic distribution of the element Al
    """
    distribution_path = file_path("/Al/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Al"].value)


def test_electronic_distribution_si(file_path):
    """
    Test eletronic distribution of the element Si
    """
    distribution_path = file_path("/Si/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Si"].value)


def test_electronic_distribution_p(file_path):
    """
    Test eletronic distribution of the element P
    """
    distribution_path = file_path("/P/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["P"].value)


def test_electronic_distribution_s(file_path):
    """
    Test eletronic distribution of the element S
    """
    distribution_path = file_path("/S/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["S"].value)


def test_electronic_distribution_cl(file_path):
    """
    Test eletronic distribution of the element Cl
    """
    distribution_path = file_path("/Cl/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Cl"].value)


def test_electronic_distribution_ar(file_path):
    """
    Test eletronic distribution of the element Ar
    """
    distribution_path = file_path("/Ar/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Ar"].value)


def test_electronic_distribution_k(file_path):
    """
    Test eletronic distribution of the element K
    """
    distribution_path = file_path("/K/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["K"].value)


def test_electronic_distribution_ca(file_path):
    """
    Test eletronic distribution of the element Ca
    """
    distribution_path = file_path("/Ca/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Ca"].value)


def test_electronic_distribution_sc(file_path):
    """
    Test eletronic distribution of the element Sc
    """
    distribution_path = file_path("/Sc/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Sc"].value)


def test_electronic_distribution_ti(file_path):
    """
    Test eletronic distribution of the element Ti
    """
    distribution_path = file_path("/Ti/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Ti"].value)


def test_electronic_distribution_v(file_path):
    """
    Test eletronic distribution of the element V
    """
    distribution_path = file_path("/V/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["V"].value)


def test_electronic_distribution_cr(file_path):
    """
    Test eletronic distribution of the element Cr
    """
    distribution_path = file_path("/Cr/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Cr"].value)


def test_electronic_distribution_mn(file_path):
    """
    Test eletronic distribution of the element Mn
    """
    distribution_path = file_path("/Mn/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Mn"].value)


def test_electronic_distribution_fe(file_path):
    """
    Test eletronic distribution of the element Fe
    """
    distribution_path = file_path("/Fe/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Fe"].value)


def test_electronic_distribution_co(file_path):
    """
    Test eletronic distribution of the element Co
    """
    distribution_path = file_path("/Co/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Co"].value)


def test_electronic_distribution_ni(file_path):
    """
    Test eletronic distribution of the element Ni
    """
    distribution_path = file_path("/Ni/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Ni"].value)


def test_electronic_distribution_cu(file_path):
    """
    Test eletronic distribution of the element Cu
    """
    distribution_path = file_path("/Cu/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Cu"].value)


def test_electronic_distribution_zn(file_path):
    """
    Test eletronic distribution of the element Zn
    """
    distribution_path = file_path("/Zn/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Zn"].value)


def test_electronic_distribution_ga(file_path):
    """
    Test eletronic distribution of the element Ga
    """
    distribution_path = file_path("/Ga/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Ga"].value)


def test_electronic_distribution_ge(file_path):
    """
    Test eletronic distribution of the element Ge
    """
    distribution_path = file_path("/Ge/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Ge"].value)


def test_electronic_distribution_as(file_path):
    """
    Test eletronic distribution of the element As
    """
    distribution_path = file_path("/As/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["As"].value)


def test_electronic_distribution_se(file_path):
    """
    Test eletronic distribution of the element Se
    """
    distribution_path = file_path("/Se/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Se"].value)


def test_electronic_distribution_br(file_path):
    """
    Test eletronic distribution of the element Br
    """
    distribution_path = file_path("/Br/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Br"].value)


def test_electronic_distribution_kr(file_path):
    """
    Test eletronic distribution of the element Kr
    """
    distribution_path = file_path("/Kr/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Kr"].value)


def test_electronic_distribution_rb(file_path):
    """
    Test eletronic distribution of the element Rb
    """
    distribution_path = file_path("/Rb/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Rb"].value)


def test_electronic_distribution_sr(file_path):
    """
    Test eletronic distribution of the element Sr
    """
    distribution_path = file_path("/Sr/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Sr"].value)


def test_electronic_distribution_y(file_path):
    """
    Test eletronic distribution of the element Y
    """
    distribution_path = file_path("/Y/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Y"].value)


def test_electronic_distribution_zr(file_path):
    """
    Test eletronic distribution of the element Zr
    """
    distribution_path = file_path("/Zr/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Zr"].value)


def test_electronic_distribution_nb(file_path):
    """
    Test eletronic distribution of the element Nb
    """
    distribution_path = file_path("/Nb/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Nb"].value)


def test_electronic_distribution_mo(file_path):
    """
    Test eletronic distribution of the element Mo
    """
    distribution_path = file_path("/Mo/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Mo"].value)


def test_electronic_distribution_tc(file_path):
    """
    Test eletronic distribution of the element Tc
    """
    distribution_path = file_path("/Tc/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Tc"].value)


def test_electronic_distribution_ru(file_path):
    """
    Test eletronic distribution of the element Ru
    """
    distribution_path = file_path("/Ru/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Ru"].value)


def test_electronic_distribution_rh(file_path):
    """
    Test eletronic distribution of the element Rh
    """
    distribution_path = file_path("/Rh/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Rh"].value)


def test_electronic_distribution_pd(file_path):
    """
    Test eletronic distribution of the element Pd
    """
    distribution_path = file_path("/Pd/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Pd"].value)


def test_electronic_distribution_ag(file_path):
    """
    Test eletronic distribution of the element Ag
    """
    distribution_path = file_path("/Ag/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Ag"].value)


def test_electronic_distribution_cd(file_path):
    """
    Test eletronic distribution of the element Cd
    """
    distribution_path = file_path("/Cd/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Cd"].value)


def test_electronic_distribution_in(file_path):
    """
    Test eletronic distribution of the element In
    """
    distribution_path = file_path("/In/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["In"].value)


def test_electronic_distribution_sn(file_path):
    """
    Test eletronic distribution of the element Sn
    """
    distribution_path = file_path("/Sn/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Sn"].value)


def test_electronic_distribution_sb(file_path):
    """
    Test eletronic distribution of the element Sb
    """
    distribution_path = file_path("/Sb/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Sb"].value)


def test_electronic_distribution_i(file_path):
    """
    Test eletronic distribution of the element I
    """
    distribution_path = file_path("/I/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["I"].value)


def test_electronic_distribution_xe(file_path):
    """
    Test eletronic distribution of the element Xe
    """
    distribution_path = file_path("/Xe/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Xe"].value)


def test_electronic_distribution_cs(file_path):
    """
    Test eletronic distribution of the element Cs
    """
    distribution_path = file_path("/Cs/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Cs"].value)


def test_electronic_distribution_ba(file_path):
    """
    Test eletronic distribution of the element Ba
    """
    distribution_path = file_path("/Ba/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Ba"].value)


def test_electronic_distribution_ce(file_path):
    """
    Test eletronic distribution of the element Ce
    """
    distribution_path = file_path("/Ce/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Ce"].value)


def test_electronic_distribution_nd(file_path):
    """
    Test eletronic distribution of the element Nd
    """
    distribution_path = file_path("/Nd/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Nd"].value)


def test_electronic_distribution_pm(file_path):
    """
    Test eletronic distribution of the element Pm
    """
    distribution_path = file_path("/Pm/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Pm"].value)


def test_electronic_distribution_sm(file_path):
    """
    Test eletronic distribution of the element Sm
    """
    distribution_path = file_path("/Sm/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Sm"].value)


def test_electronic_distribution_gd(file_path):
    """
    Test eletronic distribution of the element Gd
    """
    distribution_path = file_path("/Gd/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Gd"].value)


def test_electronic_distribution_er(file_path):
    """
    Test eletronic distribution of the element Er
    """
    distribution_path = file_path("/Er/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Er"].value)


def test_electronic_distribution_tm(file_path):
    """
    Test eletronic distribution of the element Tm
    """
    distribution_path = file_path("/Tm/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Tm"].value)


def test_electronic_distribution_yb(file_path):
    """
    Test eletronic distribution of the element Yb
    """
    distribution_path = file_path("/Yb/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Yb"].value)


def test_electronic_distribution_lu(file_path):
    """
    Test eletronic distribution of the element Lu
    """
    distribution_path = file_path("/Lu/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Lu"].value)


def test_electronic_distribution_hf(file_path):
    """
    Test eletronic distribution of the element Hf
    """
    distribution_path = file_path("/Hf/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Hf"].value)


def test_electronic_distribution_ta(file_path):
    """
    Test eletronic distribution of the element Ta
    """
    distribution_path = file_path("/Ta/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Ta"].value)


def test_electronic_distribution_w(file_path):
    """
    Test eletronic distribution of the element W
    """
    distribution_path = file_path("/W/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["W"].value)


def test_electronic_distribution_re(file_path):
    """
    Test eletronic distribution of the element Re
    """
    distribution_path = file_path("/Re/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Re"].value)


def test_electronic_distribution_os(file_path):
    """
    Test eletronic distribution of the element Os
    """
    distribution_path = file_path("/Os/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Os"].value)


def test_electronic_distribution_ir(file_path):
    """
    Test eletronic distribution of the element Ir
    """
    distribution_path = file_path("/Ir/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Ir"].value)


def test_electronic_distribution_pt(file_path):
    """
    Test eletronic distribution of the element Pt
    """
    distribution_path = file_path("/Pt/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Pt"].value)


def test_electronic_distribution_au(file_path):
    """
    Test eletronic distribution of the element Au
    """
    distribution_path = file_path("/Au/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Au"].value)


def test_electronic_distribution_hg(file_path):
    """
    Test eletronic distribution of the element Hg
    """
    distribution_path = file_path("/Hg/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Hg"].value)


def test_electronic_distribution_tl(file_path):
    """
    Test eletronic distribution of the element Tl
    """
    distribution_path = file_path("/Tl/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Tl"].value)


def test_electronic_distribution_pb(file_path):
    """
    Test eletronic distribution of the element Pb
    """
    distribution_path = file_path("/Pb/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Pb"].value)


def test_electronic_distribution_bi(file_path):
    """
    Test eletronic distribution of the element Bi
    """
    distribution_path = file_path("/Bi/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Bi"].value)


def test_electronic_distribution_po(file_path):
    """
    Test eletronic distribution of the element Po
    """
    distribution_path = file_path("/Po/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Po"].value)


def test_electronic_distribution_at(file_path):
    """
    Test eletronic distribution of the element At
    """
    distribution_path = file_path("/At/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["At"].value)


def test_electronic_distribution_rn(file_path):
    """
    Test eletronic distribution of the element Rn
    """
    distribution_path = file_path("/Rn/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Rn"].value)


def test_electronic_distribution_ra(file_path):
    """
    Test eletronic distribution of the element Ra
    """
    distribution_path = file_path("/Ra/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Ra"].value)


def test_electronic_distribution_ac(file_path):
    """
    Test eletronic distribution of the element Ac
    """
    distribution_path = file_path("/Ac/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Ac"].value)


def test_electronic_distribution_th(file_path):
    """
    Test eletronic distribution of the element Th
    """
    distribution_path = file_path("Th/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Th"].value)


def test_electronic_distribution_u(file_path):
    """
    Test eletronic distribution of the element U
    """
    distribution_path = file_path("U/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["U"].value)


def test_electronic_distribution_np(file_path):
    """
    Test eletronic distribution of the element Np
    """
    distribution_path = file_path("Np/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Np"].value)


def test_electronic_distribution_pu(file_path):
    """
    Test eletronic distribution of the element Pu
    """
    distribution_path = file_path("Pu/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Pu"].value)


def test_electronic_distribution_am(file_path):
    """
    Test eletronic distribution of the element Am
    """
    distribution_path = file_path("Am/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Am"].value)


def test_electronic_distribution_cm(file_path):
    """
    Test eletronic distribution of the element Cm
    """
    distribution_path = file_path("Cm/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Cm"].value)


def test_electronic_distribution_bk(file_path):
    """
    Test eletronic distribution of the element Bk
    """
    distribution_path = file_path("Bk/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Bk"].value)


def test_electronic_distribution_cf(file_path):
    """
    Test eletronic distribution of the element Cf
    """
    distribution_path = file_path("Cf/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Cf"].value)


def test_electronic_distribution_es(file_path):
    """
    Test eletronic distribution of the element Es
    """
    distribution_path = file_path("Es/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Es"].value)


def test_electronic_distribution_fm(file_path):
    """
    Test eletronic distribution of the element Fm
    """
    distribution_path = file_path("Fm/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Fm"].value)


def test_electronic_distribution_md(file_path):
    """
    Test eletronic distribution of the element Md
    """
    distribution_path = file_path("Md/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Md"].value)


def test_electronic_distribution_no(file_path):
    """
    Test eletronic distribution of the element No
    """
    distribution_path = file_path("No/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["No"].value)


def test_electronic_distribution_lr(file_path):
    """
    Test eletronic distribution of the element Lr
    """
    distribution_path = file_path("Lr/electronic_distribution.txt")

    with open(distribution_path, "r") as file:
        lines = file.readlines()
        assert Counter(lines) == Counter(ElectronicDistribution["Lr"].value)
