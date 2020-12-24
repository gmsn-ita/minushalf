"""
List exchange correlations code
for INP file
"""
from enum import Enum, unique


@unique
class ExchangeCorreltion(Enum):
    """
    Enum type for the exchange correlations code
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

    def __str__(self):
        return str(self.name)

    @staticmethod
    def to_list():
        """
        Generate list of exchange correlation types
        """
        return list(map(lambda element: element.value, ExchangeCorreltion))
