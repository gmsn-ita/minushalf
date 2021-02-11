"""
Enum type for correction code used
in minushalf.yaml
"""
from enum import Enum, unique


@unique
class CorrectionCode(Enum):
    """
    Enum type for the softwares suport by the program
    """

    v = "v"  # Simple valence correction
    vf = "vf"  # Fractionary valence correction
    vc = "vc"  # Simple valence and simple conduction corrections
    vfc = "vfc"  # Valence fractionary and simple conduction corrections
    vcf = "vcf"  # Simple valence  and conduction fractionary corrections
    vfcf = "vfcf"  # Valence fractionary and conduction fractionary corrections

    def __str__(self):
        return str(self.name)

    @staticmethod
    def to_list():
        """
        Generate list of available softwares
        """
        return list(map(lambda element: element.value, CorrectionCode))
