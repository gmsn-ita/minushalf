"""
List exchange and correlation codes for the INP file
"""
from enum import Enum, unique


@unique
class ExchangeCorrelation(Enum):
    """
    Enum type for exchange and correlation codes
    """

    ca = "ca"  # Ceperley-Alder
    wi = "wi"  # Wigner
    hl = "hl"  # Hedin-Lundqvist
    gl = "gl"  # Gunnarson-Lundqvist
    bh = "bh"  # Von Barth-Hedin
    pb = "pb"  # PBE scheme by Perdew, Burke, and Ernzerhof
    rp = "rp"  # RPBE scheme by Hammer, Hansen, and Norskov
    rv = "rv"  # revPBE scheme by Zhang and Yang
    bl = "bl"  # BLYP (Becke-Lee-Yang-Parr) scheme
    pz = "pz"  # Perdew-Zunger LDA

    def __str__(self):
        return str(self.name)

    @staticmethod
    def get_default():
        """
        Returns the default value for this parameter
        """
        return ExchangeCorrelation.pb.value

    @staticmethod
    def to_list():
        """
        Generate list of exchange and correlation codes
        """
        return list(map(lambda element: element.value, ExchangeCorrelation))

@unique
class ExchangeCorrelationQE(Enum):
    """
    Maps exchange and correlation codes to their QE ld1.x equivalents.
    Keys match ExchangeCorrelation so lookup is consistent.
    ca and bh are VASP-only and have no QE equivalent.
    """
    wi = "wig"
    hl = "hl"
    gl = "gl"
    pb = "pbe"
    rp = "rpbe"
    rv = "revpbe"
    bl = "blyp"
    pz = "pz"
