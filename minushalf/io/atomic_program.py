"""
Class for atomic program input parameters in minushalf.yaml
"""
from minushalf.utils.calculation_code import CalculationCode
from minushalf.utils.exchange_correlation import ExchangeCorrelation
from minushalf.io.minushalf_yaml_tags_interface import MinushalfYamlTags


class AtomicProgram(MinushalfYamlTags):
    """
    Set parameters and their default values
    """
    def __init__(
        self,
        exchange_correlation_code: str = ExchangeCorrelation.get_default(),
        calculation_code: str = CalculationCode.get_default(),
        max_iterations: int = 100,
    ) -> None:
        """
            Args:
                exchange_correlation_code (str): Code for the exchange and correlation functional
                calculation_code (str):Calculation code for inp file (ae)
                max_iterations (int): Maximum number of iterations for atomic program
        """
        self.exchange_correlation_code = exchange_correlation_code
        self.calculation_code = calculation_code
        self.max_iterations = max_iterations

    def to_list(self):
        """
        return list with the class variables
        """
        return list(self.__dict__.values())

    def to_dict(self):
        """
        Return dictionary with the class variables
        """
        return self.__dict__
