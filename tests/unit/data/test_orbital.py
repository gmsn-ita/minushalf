"""
Test Obital and OrbitalType classes
"""
from minushalf.data.orbital import Orbital, OrbitalType


def test_orbital():
    """
    Test orbital class
    """
    assert Orbital.s.value == 0
    assert Orbital.py.value == 1
    assert Orbital.pz.value == 2
    assert Orbital.px.value == 3
    assert Orbital.dxy.value == 4
    assert Orbital.dyz.value == 5
    assert Orbital.dz2.value == 6
    assert Orbital.dxz.value == 7
    assert Orbital.dx2.value == 8
    assert Orbital.f_3.value == 9
    assert Orbital.f_2.value == 10
    assert Orbital.f_1.value == 11
    assert Orbital.f0.value == 12
    assert Orbital.f1.value == 13
    assert Orbital.f2.value == 14
    assert Orbital.f3.value == 15
    assert str(Orbital.s) == "s"
    assert str(Orbital.py) == "py"
    assert str(Orbital.pz) == "pz"
    assert str(Orbital.px) == "px"
    assert str(Orbital.dxy) == "dxy"
    assert str(Orbital.dyz) == "dyz"
    assert str(Orbital.dz2) == "dz2"
    assert str(Orbital.dxz) == "dxz"
    assert str(Orbital.dx2) == "dx2"
    assert str(Orbital.f_3) == "f_3"
    assert str(Orbital.f_2) == "f_2"
    assert str(Orbital.f_1) == "f_1"
    assert str(Orbital.f0) == "f0"
    assert str(Orbital.f1) == "f1"
    assert str(Orbital.f2) == "f2"
    assert str(Orbital.f3) == "f3"


def test_orbital_type():
    """
    Test orbital type class
    """
    assert OrbitalType.s.value == 0
    assert OrbitalType.p.value == 1
    assert OrbitalType.d.value == 2
    assert OrbitalType.f.value == 3
    assert str(OrbitalType.s) == "s"
    assert str(OrbitalType.p) == "p"
    assert str(OrbitalType.d) == "d"
    assert str(OrbitalType.f) == "f"
