"""
Top-level `openff-benchmark` cli entrypoint.

"""

import click

from openff.benchmark.utils.validate_and_assign_ids import validate_and_assign
from openff.benchmark.utils.generate_conformers import generate_conformers

@click.group()
def cli():
    pass


@cli.group()
def report():
    """Analyze the results and create plots

    """
    pass

@report.command()
@click.option('--input-path', default='./', multiple=True, required=True)
@click.option('--ref_method', default='default', required=True)
@click.option('--output_directory', default='./results', required=True)
def compare_forcefields(input_path, ref_method, output_directory):
    from .analysis import analysis
    analysis.main(input_path, ref_method, output_directory)



@report.command()
@click.option('--input-path', default='./', multiple=True, required=True)
def plots(input_path):
    from .analysis import draw
    draw.plot_compare_ffs(input_path)
    if input_graph_molecules == '':
        input_graph_molecules = []
    else:
        input_graph_molecules = [input_graph_molecules]

    if input_3d_molecules == '':
        input_3d_molecules = []
    else:
        input_3d_molecules = [input_3d_molecules]

        validate_and_assign(input_graph_molecules,
                        input_3d_molecules,
                        group_name,
                        output_directory)

@cli.group()
def generate_conformers():
    pass



@generate_conformers.command()

@click.option('--input_directory',
              default='',
              help='Output directory from the validate step')

@click.option('--output_directory',
              default='2-generate_conformers', 
              help='Directory to put output files. If this directory does not exist, one will be created.')

@click.option('--group_name',
              help='Group name for assigning IDs to the molecules.')
def doit(input_directory,
         output_directory):
    generate_conformers(input_directory,
                        output_directory)
    
if __name__ == "__main__":
    cli()
