"""
Default cli messages
"""
import click
from pyfiglet import Figlet


def welcome_message(text: str) -> None:
    """
    Print welcome message
    """
    renderizer = Figlet()
    click.echo(renderizer.renderText(text))


def end_message() -> None:
    """
    Print end message
    """
    renderizer = Figlet()
    click.echo(renderizer.renderText("END"))
