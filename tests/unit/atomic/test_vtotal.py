"""
Test Vtotal class
"""
import numpy as np
from minushalf.atomic import Vtotal


def test_vtotal_cl(file_path):
    """
    Test VTOTAL parser for Cloro
    """
    vtotal_path = file_path("/Cl/VTOTAL")
    vtotal_cl = Vtotal.from_file(vtotal_path)

    assert len(vtotal_cl.radius) == len(vtotal_cl.down_potential)
    assert np.isclose(0.1834e-5, vtotal_cl.radius[0])
    assert np.isclose(119.019272152, vtotal_cl.radius[-1])
    assert np.isclose(-33.9680314505, vtotal_cl.down_potential[0])
    assert np.isclose(-0.422869e-01, vtotal_cl.down_potential[-1])


def test_vtotal_zn(file_path):
    """
    Test VTOTAL parser for Zinco
    """
    vtotal_path = file_path("/Zn/VTOTAL")
    vtotal_zn = Vtotal.from_file(vtotal_path)

    assert len(vtotal_zn.radius) == len(vtotal_zn.down_potential)
    assert np.isclose(0.10392e-5, vtotal_zn.radius[0])
    assert np.isclose(119.857307249, vtotal_zn.radius[-1])
    assert np.isclose(-59.9670386149, vtotal_zn.down_potential[0])
    assert np.isclose(-0.425847e-01, vtotal_zn.down_potential[-1])


def test_vtotal_ti(file_path):
    """
    Test VTOTAL parser for Titanium
    """
    vtotal_path = file_path("/Ti/VTOTAL")
    vtotal_ti = Vtotal.from_file(vtotal_path)

    assert len(vtotal_ti.radius) == len(vtotal_ti.down_potential)
    assert np.isclose(0.141722e-05, vtotal_ti.radius[0])
    assert np.isclose(119.576532343, vtotal_ti.radius[-1])
    assert np.isclose(-43.9675322335, vtotal_ti.down_potential[0])
    assert np.isclose(-0.4248498e-01, vtotal_ti.down_potential[-1])


def test_vtotal_p(file_path):
    """
    Test VTOTAL parser for Phosphor
    """
    vtotal_path = file_path("/P/VTOTAL")
    vtotal_p = Vtotal.from_file(vtotal_path)

    assert len(vtotal_p.radius) == len(vtotal_p.down_potential)
    assert np.isclose(0.20785e-5, vtotal_p.radius[0])
    assert np.isclose(119.038671474, vtotal_p.radius[-1])
    assert np.isclose(-29.9683278539, vtotal_p.down_potential[0])
    assert np.isclose(-0.4229388e-01, vtotal_p.down_potential[-1])


def test_vtotal_c(file_path):
    """
    Test VTOTAL parser for Carbon
    """
    vtotal_path = file_path("/C/VTOTAL")
    vtotal_c = Vtotal.from_file(vtotal_path)

    assert len(vtotal_c.radius) == len(vtotal_c.down_potential)
    assert np.isclose(0.5196477e-5, vtotal_c.radius[0])
    assert np.isclose(119.490524267, vtotal_c.radius[-1])
    assert np.isclose(-11.9720173404, vtotal_c.down_potential[0])
    assert np.isclose(-0.4245442e-01, vtotal_c.down_potential[-1])


def test_vtotal_bi(file_path):
    """
    Test VTOTAL parser for Bismuth
    """
    vtotal_path = file_path("/Bi/VTOTAL")
    vtotal_bi = Vtotal.from_file(vtotal_path)

    assert len(vtotal_bi.radius) == len(vtotal_bi.down_potential)
    assert np.isclose(0.375648e-6, vtotal_bi.radius[0])
    assert np.isclose(119.242488916, vtotal_bi.radius[-1])
    assert np.isclose(-165.966198780, vtotal_bi.down_potential[0])
    assert np.isclose(-0.42366300e-01, vtotal_bi.down_potential[-1])


def test_vtotal_ni(file_path):
    """
    Test VTOTAL parser for Nickel
    """
    vtotal_path = file_path("/Ni/VTOTAL")
    vtotal_ni = Vtotal.from_file(vtotal_path)

    assert len(vtotal_ni.radius) == len(vtotal_ni.down_potential)
    assert np.isclose(0.111353e-5, vtotal_ni.radius[0])
    assert np.isclose(119.139460842, vtotal_ni.radius[-1])
    assert np.isclose(-55.9671391768, vtotal_ni.down_potential[0])
    assert np.isclose(-0.42329694e-1, vtotal_ni.down_potential[-1])


def test_vtotal_k(file_path):
    """
    Test VTOTAL parser for Potassium
    """
    vtotal_path = file_path("/K/VTOTAL")
    vtotal_k = Vtotal.from_file(vtotal_path)

    assert len(vtotal_k.radius) == len(vtotal_k.down_potential)
    assert np.isclose(0.16409928e-5, vtotal_k.radius[0])
    assert np.isclose(119.171058353, vtotal_k.radius[-1])
    assert np.isclose(-37.9677993761, vtotal_k.down_potential[0])
    assert np.isclose(-0.423409214e-1, vtotal_k.down_potential[-1])


def test_vtotal_i(file_path):
    """
    Test VTOTAL parser for Iodine
    """
    vtotal_path = file_path("/I/VTOTAL")
    vtotal_i = Vtotal.from_file(vtotal_path)

    assert len(vtotal_i.radius) == len(vtotal_i.down_potential)
    assert np.isclose(0.5882804e-6, vtotal_i.radius[0])
    assert np.isclose(119.069540176, vtotal_i.radius[-1])
    assert np.isclose(-105.966477381, vtotal_i.down_potential[0])
    assert np.isclose(-0.4230485e-1, vtotal_i.down_potential[-1])


def test_vtotal_po(file_path):
    """
    Test VTOTAL parser for Polonium
    """
    vtotal_path = file_path("/Po/VTOTAL")
    vtotal_po = Vtotal.from_file(vtotal_path)

    assert len(vtotal_po.radius) == len(vtotal_po.down_potential)
    assert np.isclose(0.3711769e-6, vtotal_po.radius[0])
    assert np.isclose(119.304965932, vtotal_po.radius[-1])
    assert np.isclose(-167.966193113, vtotal_po.down_potential[0])
    assert np.isclose(-0.4238849e-1, vtotal_po.down_potential[-1])


def test_vtotal_mn(file_path):
    """
    Test VTOTAL parser for Manganese
    """
    vtotal_path = file_path("/Mn/VTOTAL")
    vtotal_mn = Vtotal.from_file(vtotal_path)

    assert len(vtotal_mn.radius) == len(vtotal_mn.down_potential)
    assert np.isclose(0.12471545e-5, vtotal_mn.radius[0])
    assert np.isclose(119.238220332, vtotal_mn.radius[-1])
    assert np.isclose(-49.9673154957, vtotal_mn.down_potential[0])
    assert np.isclose(-0.423647837e-1, vtotal_mn.down_potential[-1])


def test_vtotal_n(file_path):
    """
    Test VTOTAL parser for Nitrogen
    """
    vtotal_path = file_path("/N/VTOTAL")
    vtotal_n = Vtotal.from_file(vtotal_path)

    assert len(vtotal_n.radius) == len(vtotal_n.down_potential)
    assert np.isclose(0.445412344e-05, vtotal_n.radius[0])
    assert np.isclose(118.995642542, vtotal_n.radius[-1])
    assert np.isclose(-13.9711782037, vtotal_n.down_potential[0])
    assert np.isclose(-0.42278597e-1, vtotal_n.down_potential[-1])


def test_vtotal_h(file_path):
    """
    Test VTOTAL parser for Hydrogen
    """
    vtotal_path = file_path("/H/VTOTAL")
    vtotal_h = Vtotal.from_file(vtotal_path)

    assert len(vtotal_h.radius) == len(vtotal_h.down_potential)
    assert np.isclose(0.31178864e-4, vtotal_h.radius[0])
    assert np.isclose(119.998512118, vtotal_h.radius[-1])
    assert np.isclose(-1.99796354987, vtotal_h.down_potential[0])
    assert np.isclose(-0.4263491e-1, vtotal_h.down_potential[-1])


def test_vtotal_cu(file_path):
    """
    Test VTOTAL parser for Cooper
    """
    vtotal_path = file_path("/Cu/VTOTAL")
    vtotal_cu = Vtotal.from_file(vtotal_path)

    assert len(vtotal_cu.radius) == len(vtotal_cu.down_potential)
    assert np.isclose(0.107513324e-5, vtotal_cu.radius[0])
    assert np.isclose(119.426778853, vtotal_cu.radius[-1])
    assert np.isclose(-57.9670822880, vtotal_cu.down_potential[0])
    assert np.isclose(-0.42431777e-1, vtotal_cu.down_potential[-1])


def test_vtotal_cs(file_path):
    """
    Test VTOTAL parser for Cesium
    """
    vtotal_path = file_path("/Cs/VTOTAL")
    vtotal_cs = Vtotal.from_file(vtotal_path)

    assert len(vtotal_cs.radius) == len(vtotal_cs.down_potential)
    assert np.isclose(0.5668884e-6, vtotal_cs.radius[0])
    assert np.isclose(119.124174998, vtotal_cs.radius[-1])
    assert np.isclose(-109.966452133, vtotal_cs.down_potential[0])
    assert np.isclose(-0.423242639e-1, vtotal_cs.down_potential[-1])


def test_vtotal_sb(file_path):
    """
    Test VTOTAL parser for Antimony
    """
    vtotal_path = file_path("/Sb/VTOTAL")
    vtotal_sb = Vtotal.from_file(vtotal_path)

    assert len(vtotal_sb.radius) == len(vtotal_sb.down_potential)
    assert np.isclose(0.61135027e-6, vtotal_sb.radius[0])
    assert np.isclose(119.184648607, vtotal_sb.radius[-1])
    assert np.isclose(-101.966505327, vtotal_sb.down_potential[0])
    assert np.isclose(-0.42345749e-1, vtotal_sb.down_potential[-1])


def test_vtotal_s(file_path):
    """
    Test VTOTAL parser for Sulfur
    """
    vtotal_path = file_path("/S/VTOTAL")
    vtotal_s = Vtotal.from_file(vtotal_path)

    assert len(vtotal_s.radius) == len(vtotal_s.down_potential)
    assert np.isclose(0.19486790e-5, vtotal_s.radius[0])
    assert np.isclose(118.796265787, vtotal_s.radius[-1])
    assert np.isclose(-31.9681713318, vtotal_s.down_potential[0])
    assert np.isclose(-0.42207759e-1, vtotal_s.down_potential[-1])


def test_vtotal_tm(file_path):
    """
    Test VTOTAL parser for Thulium
    """
    vtotal_path = file_path("/Tm/VTOTAL")
    vtotal_tm = Vtotal.from_file(vtotal_path)

    assert len(vtotal_tm.radius) == len(vtotal_tm.down_potential)
    assert np.isclose(0.45186759e-6, vtotal_tm.radius[0])
    assert np.isclose(118.913126048, vtotal_tm.radius[-1])
    assert np.isclose(-137.966300055, vtotal_tm.down_potential[0])
    assert np.isclose(-0.4224927929e-1, vtotal_tm.down_potential[-1])


def test_vtotal_rb(file_path):
    """
    Test VTOTAL parser for Rubidium
    """
    vtotal_path = file_path("/Rb/VTOTAL")
    vtotal_rb = Vtotal.from_file(vtotal_path)

    assert len(vtotal_rb.radius) == len(vtotal_rb.down_potential)
    assert np.isclose(0.842672003e-6, vtotal_rb.radius[0])
    assert np.isclose(118.697889692, vtotal_rb.radius[-1])
    assert np.isclose(-73.9667917002, vtotal_rb.down_potential[0])
    assert np.isclose(-0.42172806e-1, vtotal_rb.down_potential[-1])


def test_vtotal_pb(file_path):
    """
    Test VTOTAL parser for Lead
    """
    vtotal_path = file_path("/Pb/VTOTAL")
    vtotal_pb = Vtotal.from_file(vtotal_path)

    assert len(vtotal_pb.radius) == len(vtotal_pb.down_potential)
    assert np.isclose(0.38023005e-6, vtotal_pb.radius[0])
    assert np.isclose(119.197347174, vtotal_pb.radius[-1])
    assert np.isclose(-163.966204585, vtotal_pb.down_potential[0])
    assert np.isclose(-0.423502617e-1, vtotal_pb.down_potential[-1])
