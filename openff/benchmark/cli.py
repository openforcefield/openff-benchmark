"""
Top-level `openff-benchmark` cli entrypoint.

"""

import click

@click.group()
def cli():
    pass

@cli.group()
def analyze():
    """Analyze the results after energy minimizations

    """
    pass

@analyze.command()
@click.option('--input-path', default='./', multiple=True, required=True)
@click.option('--ref_method', default='default', required=True)
@click.option('--output_directory', default='./results', required=True)
def compare_forcefields(input_path, ref_method, output_directory):
    from .analysis import analysis
    analysis.main(input_path, ref_method, output_directory)


@cli.group()
def report():
    """Create plots with results. 

    """
    pass

@analyze.command()
@click.option('--input-path', default='./', required=True)
@click.option('--dataset-name', required=True)
@click.option('--season', required=True)
def plots(season, dataset_name):
    from .analysis import plotting
    plotting.main(
            server_uri, input_path, season, dataset_name)

