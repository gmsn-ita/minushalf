"""
Handler to get software object
"""
import loguru
from minushalf.softwares import Vasp
from minushalf.data import Softwares
from minushalf.handlers import BaseHandler


class CreateSoftwareModule(BaseHandler):
    """
    Create the object that have the informations about first principle calculations
    """
    def action(self, request: any) -> any:
        """
        Reads software name in input file and get the module to
        read the first principle calculations results
        """
        loguru.logger.info(
            "Creating module to read the first principle results")

        name_to_object = {Softwares.vasp.value: Vasp()}

        minushalf_yaml = request["minushalf_yaml"]
        name = minushalf_yaml.get_software_name()
        request["software_module"] = name_to_object[name]

        return request
