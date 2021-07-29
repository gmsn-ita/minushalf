"""
Handler that runs first principle calculations
"""
import loguru
from minushalf.handlers import BaseHandler


class RunCalculations(BaseHandler):
    """
    Run first principle calculations
    """
    def action(self, request: any) -> any:
        """
        Reads software name in input file and get the module to
        read the first principle calculations results
        """
        loguru.logger.info(
            "Creating module to read the first principle results")

        softwate_module = request["software_module"]
        minushalf_yaml = request["minushalf_yaml"]
        command = minushalf_yaml.get_command()

        runner = softwate_module.get_runner(command)

        runner.run()

        return request
