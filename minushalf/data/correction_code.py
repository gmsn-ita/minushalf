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

    v = "v"  # Only valence correction
    vf = "vf"  # fractionary valence correction
    vc = "vc"  # valence and conduction corrections
    vfc = "vfc"  # valence fractionary and conduction corrections
    vcf = "vcf"  # valence correction and conduction fractionary corrections
    vfcf = "vfcf"  # valence fractionary and confduction fractionary corrections

    def __str__(self):
        return str(self.name)

    @staticmethod
    def to_list():
        """
        Generate list of available softwares
        """
        return list(map(lambda element: element.value, CorrectionCode))
