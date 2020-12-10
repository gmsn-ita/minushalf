"""
Reads vasprun.xml file, an output of
VASP software
"""
import re
import xml.etree.ElementTree as ET


class Vasprun():
    """
	Reads vasprun.xml and store
	useful informations
	"""
    def __init__(self, filename: str):
        """
			Args:
				filename (str): name of the vasprun file in VASP
			Members:
				num_kpoints: Number of kpoints used in the simulation
				num_bands: Number of bands used in the simulation
				size_procar_header: Number of lines of the procar header in PROCAR file
		"""
        self.filename = filename
        self.fermi_energy = self._get_fermi_energy()
        self.atoms_map = self._get_atoms()

    def _get_fermi_energy(self) -> float:
        """
		Extract fermi energy from vasprun.xml file
			Returns:
				fermi_energy (float): Fermi energy
		"""
        efermi_xml = self._catch_xml_tag(tag="i", name="efermi")
        xml_tree = ET.fromstringlist(efermi_xml)
        try:
            fermi_energy = float(xml_tree.text)
        except ValueError as invalid_conversion:
            raise Exception('Vasprun parser does not find fermi energy'
                            ) from invalid_conversion
        finally:
            xml_tree.clear()
        return fermi_energy

    def _get_atoms(self) -> dict:
        """
		Extract informations about the
		atoms from vasprun.xml
			Returns:
				A dict containing the index of the atoms
				and their symbols. dict[index] =  symbol
		"""
        atoms_xml = self._catch_xml_tag(tag="array", name="atoms")

        atoms_tree = ET.fromstringlist(atoms_xml)

        try:
            atoms = {}
            for index, value in enumerate(
                    atoms_tree.find("set").findall("rc")):
                atom_index = list(value.findall("c"))
                atoms[str(index + 1)] = atom_index[0].text.strip()
        except:
            raise Exception("Vasprun parser do not found atoms informations")
        finally:
            atoms_tree.clear()
        return atoms

    def _catch_xml_tag(self, tag: str, name: str = None) -> list:
        """
		Catch element by the tag in a large xml file

		Args:
			tag (str): Tag to the element in xml

		Returns:
			captured_records (np.array): lines between the
			start tag and end tag
		"""

        start_tag_regex = re.compile(rf'^.*<{tag}\s*(name="{name}")?\s*>')
        end_tag_regex = re.compile(rf"^.*</{tag}\s*>")
        start_tag_identified = False
        xml_text = []
        with open(self.filename, "r") as vasprun:
            for line in vasprun:
                if start_tag_identified and end_tag_regex.match(line):
                    end_tag = end_tag_regex.search(line).group(0)
                    end_line_splited = line.partition(end_tag)
                    end_line = "".join(end_line_splited[:2])
                    start_tag_identified = False
                    xml_text.append(end_line)
                    break

                elif start_tag_identified:
                    xml_text.append(line)

                elif not start_tag_identified and start_tag_regex.match(line):
                    start_tag = start_tag_regex.search(line).group(0)
                    start_line_splited = line.partition(start_tag)
                    start_line = "".join(start_line_splited[1:])
                    start_tag_identified = True
                    xml_text.append(start_line)
                    if end_tag_regex.match(line):
                        break
        return xml_text
