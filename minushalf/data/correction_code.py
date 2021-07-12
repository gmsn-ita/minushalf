"""
Enum type for correction codes used
in minushalf.yaml
"""
from enum import Enum, unique


@unique
class CorrectionCode(Enum):
    """
    Enum type for the correction codes
    """

    v = "v"       # Simple valence correction
    vf = "vf"     # Fractional valence correction
    c = "c"       # Simple conduction correction
    cf = "cf"     # Fractional conduction correction
    vc = "vc"     # Simple valence and simple conduction corrections
    vfc = "vfc"   # Valence fractional and simple conduction corrections
    vcf = "vcf"   # Simple valence  and conduction fractional corrections
    vfcf = "vfcf" # Valence fractionary and conduction fractional corrections

    def __str__(self):
        return str(self.name)

    @staticmethod
    def to_list():
        """
        Generate list of available correction codes
        """
        return list(map(lambda element: element.value, CorrectionCode))
