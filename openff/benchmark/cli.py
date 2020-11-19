
"""
Top-level `openff-benchmark` cli entrypoint.
"""

import click

from openff.benchmark.utils.validate_and_assign_ids import validate_and_assign


@click.group()
def cli():
    pass

@cli.group()
def validate():
    pass

@validate.command()
@click.option('-m',
              '--molecules',
              help='SDF file to read containing input molecules.')
@click.option('--groupname',
              help='Group name for assigning IDs to the molecules.')
def doit(molecules, groupname):
    validate_and_assign(molecules, groupname)
    
if __name__ == "__main__":
    cli()
