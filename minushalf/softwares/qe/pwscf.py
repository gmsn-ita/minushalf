"""
Reads the pw.x XML data file ($prefix.xml),
an output of Quantum ESPRESSO software
"""
import re
import xml.etree.ElementTree as ET
from collections import defaultdict


class PWSCF():
    """
    Reads a pw.x XML output file (e.g. pwscf.xml) and stores
    information parsed from it.
    """

    def __init__(self, filename: str):
        """
        Args:
            filename (str): path to the pw.x XML output file ($prefix.xml)
        Members:
            filename:    stored path
            fermi_energy: Fermi energy in Hartree atomic units
            eigenvalues:  defaultdict(list) where keys are 1-based kpoint
                          indices and values are lists of eigenvalues in
                          Hartree atomic units for each band, in order
        """
        self.filename = filename
        self.fermi_energy = self._get_fermi_energy()
        self.eigenvalues = self._get_eigenvalues()

    # ------------------------------------------------------------------ #
    #  Parsers                                                             #
    # ------------------------------------------------------------------ #

    def _get_fermi_energy(self) -> float:
        """
        Extract the Fermi energy from the QE XML data file.
        Returns:
            fermi_energy (float): Fermi energy in Hartree atomic units
        """
        fermi_xml = self._catch_xml_tag(tag="fermi_energy")
        xml_tree = ET.fromstringlist(fermi_xml)
        try:
            fermi_energy = float(xml_tree.text)
        except ValueError as invalid_conversion:
            raise Exception(
                "PWSCF parser could not parse the fermi energy value"
            ) from invalid_conversion
        finally:
            xml_tree.clear()
        return fermi_energy

    def _get_eigenvalues(self) -> defaultdict:
        """
        Extract eigenvalues for every k-point and band from the QE XML file.

        Returns:
            eigenvalues (defaultdict(list)):
                { kpoint_index: [e_band1, e_band2, ...], ... }
        """
        eigenvalues = defaultdict(list)
        kpoint = 0

        ks_open_regex   = re.compile(r"^\s*<ks_energies>")
        eig_open_regex  = re.compile(r"^\s*<eigenvalues[^>]*>")
        eig_close_regex = re.compile(r"^\s*</eigenvalues\s*>")
        data_regex      = re.compile(
            r"^\s*[-+]?[0-9]*\.?[0-9]+[EeDd][-+]?\d+"
        )

        with open(self.filename, "r") as xml_file:
            inside_ks  = False
            inside_eig = False

            for line in xml_file:
                if not inside_ks:
                    if ks_open_regex.match(line):
                        inside_ks = True
                        kpoint += 1
                    continue

                if not inside_eig:
                    if eig_open_regex.match(line):
                        inside_eig = True
                    continue

                if eig_close_regex.match(line):
                    inside_eig = False
                    inside_ks  = False
                    continue

                if data_regex.match(line):
                    eigenvalues[kpoint].extend(
                        float(v) for v in line.split()
                    )

        if not eigenvalues:
            raise Exception(
                "PWSCF parser could not find any <ks_energies> blocks "
                f"in {self.filename}"
            )

        return eigenvalues

    # ------------------------------------------------------------------ #
    #  Private helpers                                                   #
    # ------------------------------------------------------------------ #

    def _catch_xml_tag(self, tag: str, name: str = None) -> list:
        """
        Capture the lines between an opening and closing XML tag in a
        large XML file without loading the whole document into memory.

        Args:
            tag (str):  XML tag name to search for (e.g. 'fermi_energy')
            name (str): optional value of a 'name' attribute used to
                        disambiguate tags that share the same tag name
                        (mirrors the Vasprun helper signature)
        Returns:
            xml_text (list[str]): lines from the opening tag to the
                                  closing tag, inclusive
        """
        start_tag_regex = re.compile(
            rf'^.*<{tag}\s*(name="{name}")?\s*>'
        )
        end_tag_regex = re.compile(rf'^.*</{tag}\s*>')

        xml_text = []
        with open(self.filename, "r") as xml_file:
            start_tag_identified = False
            for line in xml_file:
                if start_tag_identified and end_tag_regex.match(line):
                    end_tag = end_tag_regex.search(line).group(0)
                    end_line_split = line.partition(end_tag)
                    xml_text.append("".join(end_line_split[:2]))
                    break
                elif start_tag_identified:
                    xml_text.append(line)
                elif not start_tag_identified and start_tag_regex.match(line):
                    start_tag = start_tag_regex.search(line).group(0)
                    start_line_split = line.partition(start_tag)
                    xml_text.append("".join(start_line_split[1:]))
                    start_tag_identified = True
                    if end_tag_regex.match(line):
                        break

        if not xml_text:
            raise Exception(
                f"PWSCF parser could not find tag '<{tag}>' in {self.filename}"
            )

        return xml_text
    