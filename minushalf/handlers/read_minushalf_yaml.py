"""
Handler to read minushalf.yaml
"""
from minushalf.io import MinushalfYaml
from .base_handler import BaseHandler
from loguru import logger


class ReadMinushalf(BaseHandler):
    """
    Read minushalf yaml and pass it to other handlers
    """

    def action(self, request: any) -> any:
        """
        Reads minushalf.yaml and return it
        """
        logger.info("Reading minushalf.yaml file")

        minushalf_yaml = MinushalfYaml.from_file()
        request['minushalf_yaml'] = minushalf_yaml

        return request