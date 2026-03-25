"""plantvelo.commands
~~~~~~~~~~~~~~~~~~~~
CLI entry point: ``plantvelo``
"""
import click
import logging
import sys
from collections import OrderedDict

from plantvelo.commands.run import run
import plantvelo._version


class NaturalOrderGroup(click.Group):
    """Preserve registration order in --help output."""

    def list_commands(self, ctx):
        return self.commands.keys()


@click.version_option(version=plantvelo._version.__version__)
@click.group(
    cls=NaturalOrderGroup,
    commands=OrderedDict(),
    context_settings=dict(max_content_width=300, terminal_width=300),
)
def cli() -> None:
    """PlantVelocity – RNA velocity for plant single-cell data.

    Extends velocyto with an intron_retained layer for three-state kinetic
    modelling (U → IR → S).
    """
    logging.basicConfig(
        stream=sys.stdout,
        format='%(asctime)s - %(levelname)s - %(message)s',
        level=logging.DEBUG,
    )


cli.add_command(run)
