"""
Vasp correction factory
"""
from minushalf.interfaces import CorrectionAbstractFactory
from .vasp import VaspValenceCorrection


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
