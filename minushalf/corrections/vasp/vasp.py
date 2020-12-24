"""
Correction of pseudofile
implemented for VASP software.
"""
from datetime import date
import numpy as np
from minushalf.corrections.base_software import Software


class Vasp(Software):
    """
    Read inputs from atomic program
    and correct pseudopotential file (POTCAR)
    for fractional occupation.
    """
    def __init__(self,
                 cut: float,
                 amplitude: float = 1.0,
                 potcar_path: str = "POTCAR"):
        """
        Args:
            cut (float): Value of limit radius necessary for DFT -1/2 techinique
            amplitude (float):  Triming function modulator
            potcar_path (str): Path to POTCAR file
        """
        super().__init__(cut, potcar_path, amplitude)
        self.potcar_first_line = ""
        self.psctr = []
        self.potcar_last_lines = np.array([])
        if np.isclose(abs(self.amplitude), 1.0):
            self.new_pot_name = "POTCARcut{:.2f}".format(self.cut)
        else:
            self.new_pot_name = "POTCARcut{:.2f}A{:.1f}".format(
                self.cut, self.amplitude)

    def parse_potfile(self):
        """
        Structure informations contained in POTCAR file and save
        in the following stuctures:

        potcar_first_line: stores the first line of POTCAR file

        psctr: stores the PSCTR params of POTCAR file

        k_max: maximum fourier coefficient

        fourier_coefficients: stores the fourier coeficients present in POTCAR

        last_lines: stores the aditional lines that will not be used in the calculations.

        """

        with open(self.potfile_path, "r") as potcar:
            # Parsing POTCAR
            potcar_lines = potcar.readlines()

            # First line
            self.potcar_first_line = potcar_lines[0].rstrip('\n')

            # psctr
            initial_psctr_index = 1
            final_psctr_index = 2 + find_element(
                potcar_lines, "END of PSCTR-controll parameters",
                initial_psctr_index)

            self.psctr = potcar_lines[initial_psctr_index:final_psctr_index]

            # k max
            k_max_index = final_psctr_index
            self.k_max = float(potcar_lines[k_max_index])

            # Fourier coeficientes
            intial_fourier_coef_index = final_psctr_index + 1
            final_fourier_coef_index = find_element(
                potcar_lines, "gradient corrections used for XC",
                intial_fourier_coef_index)
            fourier_coef_lines = potcar_lines[
                intial_fourier_coef_index:final_fourier_coef_index]
            self.fourier_coef = np.array(" ".join(fourier_coef_lines).split(),
                                         dtype=np.float)

            # Restant lines
            self.potcar_last_lines = potcar_lines[final_fourier_coef_index:]

    def write_new_potfile(self):
        """
        Write corrected POTCAR file.
        """
        with open(self.new_pot_name, "w") as new_potcar:
            new_potcar.write("{} CUT={:.3} {}Amp={:.1} LDA-1/2 {}\n".format(
                self.potcar_first_line, self.cut, self.resume_inp,
                self.amplitude, date.today()))

            new_potcar.writelines(self.psctr)
            new_potcar.write("{:.18}\n".format(self.k_max))

            for index, value in enumerate(self.corrected_fourier_coef):
                new_potcar.write(" {}".format(format(value, "15.8E")))

                if index > 0 and (index + 1) % 5 == 0:
                    new_potcar.write('\n')

            new_potcar.writelines(self.potcar_last_lines)
