
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
@click.option('--input_molecules',
              default='',
              help="Molecules to put through the standard benchmarking workflow. 3D conformers in this file will be ignored and only the molecule graph/connectivity will be used for the workflow.")

@click.option('--input_conformers',
              default='',
              help='SDF file to read containing input molecules in specific geometries. This argument should be included if there are particular conformers')

@click.option('--output_directory',
              default='1-validate_and_assign', 
              help='Directory to put output files. If this directory does not exist, one will be created.')

@click.option('--group_name',
              help='Group name for assigning IDs to the molecules.')
def doit(input_molecules,
         input_conformers,
         output_directory,
         group_name):
    validate_and_assign(input_molecules,
                        input_conformers,
                        output_directory,
                        group_name)
    
if __name__ == "__main__":
    cli()
