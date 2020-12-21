"""
Execute command
"""
import click
import yaml


@click.command()
def execute():
    """
        Execute command
            Requires:
                minushalf.yaml file : Parameters file
            Return:
                minushalf_results.dat : Results file
    """
