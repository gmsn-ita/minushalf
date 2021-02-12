"""
Leads with input file (INP.ae)
read by atomic program.
"""
import numpy as np
import fortranformat as ff
from minushalf.data import (ElectronicDistribution, PeriodicTable,
                            ExchangeCorreltion, CalculationCode)
from .drop_comments import drop_comments
from .parse_valence_orbital_line import parse_valence_orbitals


class InputFile:
    """
    Parses input file.
    """
    def __init__(self,
                 exchange_correlation_type: str,
                 calculation_code: str,
                 chemical_symbol: str,
                 esoteric_line: str,
                 number_valence_orbitals: int,
                 number_core_orbitals: int,
                 valence_orbitals: list,
                 description: str = "",
                 last_lines: list = None) -> None:
        """
        Args:
            chemical_symbol (str): Symbol of the chemical element (H, He, Li...)

            esoteric_line (str):  Its use is somewhat esoteric and for most
            calculations it should contain just a 0.0 in the position shown.

            exchange_correlation_type (str): functional of exchange and correlation
            ((r)ca(s), (r)wi(s), (r)hl(s), (r)gl(s) ,(r)bh(s), (r)pb(s), (r)rp(s), (r)rv(s), (r)bl(s))

            calculation code (str): Calculation code for inp file (ae)

            number_valence_orbitals (int): Number of orbitals in valence

            number_core_orbitals (int): Number of orbitals in the core

            valence_orbitals (list): list of dictionaries with the following
            properties: {"n": principal quantum number,"l":secondary quantum number,
            "occupation": occupation in the level}

            last_lines (list): any line or property that comes after
            electronic distribution
        """
        self.exchange_correlation_type = exchange_correlation_type
        self.calculation_code = calculation_code
        self.description = description
        self.chemical_symbol = chemical_symbol
        self.esoteric_line = esoteric_line
        self.number_core_orbitals = number_core_orbitals
        self.number_valence_orbitals = number_valence_orbitals
        self.valence_orbitals = valence_orbitals
        if not last_lines:
            self.last_lines = []
        else:
            self.last_lines = last_lines

    @property
    def chemical_symbol(self) -> str:
        """
        Returns:
            Chemical symbol of the element (H, He, Li...)
        """
        return self._chemical_symbol

    @chemical_symbol.setter
    def chemical_symbol(self, symbol: str) -> None:
        """
        Verify if the symbol is a valid periodic table element and
        format the string correctly.

        Args:
            symbol (str): chemical symbol of the element (H, He, Li...)
        """

        try:
            PeriodicTable[symbol]
        except KeyError as symbol_not_found:
            raise ValueError("The chemical symbol passed is not correct"
                             ) from symbol_not_found

        self._chemical_symbol = symbol.capitalize()

    @property
    def exchange_correlation_type(self) -> str:
        """
        Returns:
            Functional of exchange and correlation
            (ca, wi, hl, gl,bh, pb, rp, rv, bl
        """
        return self._exchange_correlation_type

    @exchange_correlation_type.setter
    def exchange_correlation_type(self,
                                  exchange_correlation_type: str) -> None:
        """
        Verify if the functional of  exchange and correlation is valid
        conforms the ATOM documentation

        Args:
            exchange_correlation_type (str): functional of exchange and correlation
            (ca, wi, hl, gl ,bh, pb, rp, rv, bl()
        """
        try:
            value = ExchangeCorreltion[exchange_correlation_type].value
        except KeyError as code_not_found:
            raise KeyError(
                "Your value of exchange and correlation functional is not valid"
            ) from code_not_found

        self._exchange_correlation_type = value

    @property
    def calculation_code(self) -> str:
        """
        Returns:
            Calculation code for inp file (ae)
        """
        return self._calculation_code

    @calculation_code.setter
    def calculation_code(self, calculation_code: str) -> None:
        """
        Verify if the calculation is valid
        conforms the ATOM documentation

        Args:
            calculation code (str): Calculation code for inp file (ae)
        """
        try:
            value = CalculationCode[calculation_code].value
        except KeyError as code_not_found:
            raise KeyError(
                "Your value of calculation is not valid") from code_not_found

        self._calculation_code = value

    def electron_occupation(self, electron_fraction: float,
                            secondary_quantum_number: int) -> None:
        """
        Corrects the input file of the atomic program,
        decreasing a fraction of the electron in a
        layer specified by the secondary quantum number

            Args:
                electron_fraction (float): Fraction of the electron
                that will be decreased in the INP file. Can vary between 0 and 0.5

                secondary_quantum_number (int): Specifies the layer on which
                the occupation is to be made.
        """
        for value in reversed(self.valence_orbitals):

            is_equal = np.isclose(value["occupation"][0],
                                  0.0,
                                  rtol=1e-04,
                                  atol=1e-08,
                                  equal_nan=False)

            if (not is_equal and value["l"] == secondary_quantum_number):
                value["occupation"][0] -= electron_fraction
                break
        else:
            raise Exception(
                "Trouble with occupation. Please verify the parameters passed and the INP file."
            )

    def to_stringlist(self) -> list:
        """
            Returns:
                List with the lines of the INP file.
        """
        input_lines = []

        input_lines.append("{:<3}{}{:<6}{}\n".format("", self.calculation_code,
                                                     "", self.description))

        chemical_symbol_line = "n={}".format(self.chemical_symbol)
        exchange_correlation_line = "c={}".format(
            self.exchange_correlation_type)

        second_line_formater = ff.FortranRecordWriter('1x,a4,1x,a4,2x')
        if len(self.chemical_symbol) == 1:
            second_line_formater = ff.FortranRecordWriter('1x,a3,2x,a4,2x')

        input_lines.append("{}\n".format(
            second_line_formater.write(
                [chemical_symbol_line, exchange_correlation_line])))

        input_lines.append(self.esoteric_line)

        orbital_numbers_formater = ff.FortranRecordWriter('2i5')
        orbital_numbers = orbital_numbers_formater.write(
            [self.number_core_orbitals, self.number_valence_orbitals])

        input_lines.append("{}\n".format(orbital_numbers))

        for orbital in self.valence_orbitals:
            quantum_number_formater = ff.FortranRecordWriter('2i5')
            quantum_numbers = quantum_number_formater.write(
                [orbital["n"], orbital["l"]])

            occupation_formater = ff.FortranRecordWriter('2f10.3')
            occupation = occupation_formater.write(orbital["occupation"])

            input_lines.append("{}{}\n".format(quantum_numbers, occupation))

        for line in self.last_lines:
            input_lines.append(line)

        return input_lines

    def to_file(self, filename: str = "./INP") -> None:
        """
        Write INP file
            Args:
                filename (str): name of the output file
        """

        with open(filename, "w") as input_file:
            lines = self.to_stringlist()
            input_file.writelines(lines)

    @staticmethod
    def from_file(filename: str = "./INP") -> any:
        """
        Parse INP.ae file.

            Args:
                filename: name of the INP file.
            Returns:
                input_file: instance of InputFile class.
        """
        with open(filename) as input_file:

            lines_without_comments = drop_comments(input_file.readlines())

            try:
                calculation_code = lines_without_comments[0].split()[0]
                description = " ".join(lines_without_comments[0].split()[1:])
            except ValueError as bad_inp_format:
                raise ValueError("Description or calculation code not provided"
                                 ) from bad_inp_format

            try:
                chemical_symbol = lines_without_comments[1].split()[0].split(
                    "=")[1]
                exchange_correlation_type = lines_without_comments[1].split(
                )[1].split("=")[1]
            except ValueError as bad_inp_format:
                raise ValueError(
                    "Chemical symbol or exchange correlation not provided"
                ) from bad_inp_format

            esoteric_line = lines_without_comments[2]

            try:
                number_core_orbitals = int(
                    lines_without_comments[3].split()[0])
                number_valence_orbitals = int(
                    lines_without_comments[3].split()[1])
            except ValueError as bad_inp_format:
                raise ValueError(
                    "Number of core orbitals or number of valence orbitals not provided"
                ) from bad_inp_format

            try:
                valence_orbitals = [
                    parse_valence_orbitals(lines_without_comments[i])
                    for i in range(4, 4 + number_valence_orbitals)
                ]
            except ValueError as bad_inp_format:
                raise ValueError("Valence orbitals do not provided correctly"
                                 ) from bad_inp_format

            last_lines = lines_without_comments[4 + number_valence_orbitals:]

            return InputFile(exchange_correlation_type, calculation_code,
                             chemical_symbol, esoteric_line,
                             number_valence_orbitals, number_core_orbitals,
                             valence_orbitals, description, last_lines)

    @staticmethod
    def minimum_setup(chemical_symbol: str,
                      exchange_correlation_type: str,
                      maximum_iterations: int = 100,
                      calculation_code: str = "ae") -> any:
        """
        Create INP file with minimum setup.

            Args:
            chemical_symbol (str): Symbol of the chemical element (H, He, Li...).

            exchange_correlation_type (str): functional of exchange and correlation
            ( ca, wi, hl, gl, bh, pb, rp , rv, bl)

            maximum_iterations (int): Maximum number of iterations for atomic program.
            The default is 100

            Returns:
                input_file: instance of InputFile class.
        """
        description = "{}".format(chemical_symbol)
        esoteric_line = "{:<7}0.0{:<7}0.0{:<7}0.0{:<7}0.0{:<7}0.0{:<7}0.0\n".format(
            '', '', '', '', '', '')
        last_lines = ["{} maxit\n".format(maximum_iterations)]

        try:
            electronic_distribution = ElectronicDistribution[
                chemical_symbol].value
        except ValueError as element_not_found:
            raise ValueError("This element its not available in our database"
                             ) from element_not_found
        number_core_orbitals = int(electronic_distribution[0].split()[0])
        number_valence_orbitals = int(electronic_distribution[0].split()[1])
        valence_orbitals = [
            parse_valence_orbitals(orbital)
            for orbital in electronic_distribution[1:]
        ]

        return InputFile(exchange_correlation_type, calculation_code,
                         chemical_symbol, esoteric_line,
                         number_valence_orbitals, number_core_orbitals,
                         valence_orbitals, description, last_lines)
