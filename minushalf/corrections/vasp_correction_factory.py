"""
Vasp correction factory
"""
from minushalf.interfaces import CorrectionAbstractFactory
from .vasp import (
    VaspValenceCorrection,
    VaspValenceFractionaryCorrection,
    VaspConductionCorrection,
    VaspConductionFractionaryCorrection,
)


class VaspCorrectionFactory(CorrectionAbstractFactory):
    """
    Factory to delivery correction
    algorithms for VASP software
    """
    def valence(self, **kwargs) -> VaspValenceCorrection:
        """
        Returns algorithm to valence correction in vasp
        """
        return VaspValenceCorrection(**kwargs)

    def valence_fractionary(self,
                            **kwargs) -> VaspValenceFractionaryCorrection:
        """
        Returns algorithm to valence fractionary correction in vasp
        """
        return VaspValenceFractionaryCorrection(**kwargs)

    def conduction(self, **kwargs) -> VaspConductionCorrection:
        """
        Returns algorithm to simple conduction correction in vasp
        """
        return VaspConductionCorrection(**kwargs)

    def conduction_fractionary(
            self, **kwargs) -> VaspConductionFractionaryCorrection:
        """
        Returns algorithm to conduction fractionary correction in vasp
        """
        return VaspConductionFractionaryCorrection(**kwargs)
