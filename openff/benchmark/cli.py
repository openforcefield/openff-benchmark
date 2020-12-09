"""
Top-level `openff-benchmark` cli entrypoint.

"""

import click

@click.group()
def cli():
    pass

@cli.group()
def fractal():
    """Control points for a fractal server and managers.

    """
    pass

@fractal.command()
def server_init():
    from .fractal import fractal_server_init
    fractal_server_init()

@fractal.command()
def server_start():
    from .fractal import fractal_server_start
    fractal_server_start()

@fractal.command()
def manager_init():
    pass

@fractal.command()
def manager_start():
    from .fractal import fractal_manager_start
    fractal_manager_start()

@cli.group()
def optimize():
    """Execute benchmarking geometry optimizations.

    """
    pass

@optimize.command()
@click.option('--server-uri', default="localhost:7777")
@click.option('--dataset-name', required=True)
#@click.option('--replace', is_flag=True)
@click.option('--season', required=True)
@click.argument('input-path', nargs=-1)
def submit_molecules(server_uri, input_path, season, dataset_name):
    from .geometry_optimizations import compute as optcompute
    optcompute.submit_molecules(
            server_uri, input_path, season, dataset_name)


@optimize.command()
def submit_compute(server_uri, input_path):
    pass

@optimize.command()
@click.option('--server-uri', default="localhost:7777")
@click.option('--dataset-name', required=True)
@click.option('--delete-existing', is_flag=True)
@click.argument('output-directory')
def export(server_uri, output_directory, dataset_name, delete_existing):
    from .geometry_optimizations import compute as optcompute
    optcompute.export_molecule_data(
            server_uri, output_directory, dataset_name, delete_existing)

@optimize.command()
@click.option('--server-uri', default="localhost:7777")
@click.option('--dataset-name', required=True)
def status(server_uri, dataset_name):
    from .geometry_optimizations import compute as optcompute
    optdf = optcompute.get_optimization_status(server_uri, dataset_name)
    print(optdf.applymap(lambda x: x.status.value))


@optimize.command()
@click.option('--mode',
              type=click.Choice(['singleshot', 'service']),
              default='singleshot',
              help="If 'singleshot', runs once and exits; if 'service', runs in a loop")
@click.option('--sleep-time', default=3600, help="Time in seconds to sleep when in 'service' mode")
@click.option('--only-compute', default=None, help="Comma-separated compute spec names to limit errorcycling to")
@click.argument('molids', nargs=-1)
def errorcycle(molids, mode, sleep_time, only_compute):
    pass

@optimize.command()
def execute_from_server():
    pass

@optimize.command()
@click.option('--season', required=True)
@click.option('--ncores', default=2)
@click.option('--delete-existing', is_flag=True)
@click.argument('input-path', nargs=-1)
@click.argument('output-directory')
def execute(input_path, output_directory, season, ncores, delete_existing):
    from .geometry_optimizations import compute as optcompute
    optcompute.execute_optimization_from_molecules(
            input_path, output_directory, season, ncores=ncores, delete_existing=delete_existing)


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



@cli.group()
def validate():
    pass

@validate.command()
@click.option('--input_graph_molecules',
              default='',
              help="SDF file(s) to read containing molecules to put through the standard benchmarking workflow. 3D conformers in this file will be ignored and the molecule graph/connectivity will be used to generate conformers for these molecules.")
@click.option('--input_3d_molecules',
              default='',
              nargs=1,
              help='SDF file(s) to read containing input molecules in specific geometries. This argument should be included if there are particular low-energy conformers that naive conformer generation may not find.')
@click.option('--group_name',
              help='Group name for assigning IDs to the molecules.')
@click.option('--output_directory',
              default='1-validate_and_assign', 
              help='Directory to put output files. If this directory does not exist, one will be created.')
def doit(input_graph_molecules,
         input_3d_molecules,
         output_directory,
         group_name):
    """
    This script preprocesses, validates, and creates a naming system for molecles that will be submitted to the benchmarking workflow. 
    For each unique molecule, up to ten total conformers will be generated by the benchmarking workflow.\n
    This script takes two forms of molecule inputs:\n
        1) "Graph Molecules", the "normal" input type where the geometry of the input molecule will NOT be passed to subsequent steps, and instead all conformers of the molecule will be generated by RDKit.\n
        2) "3D Molecules", where one or more specific geometries are provided, and the rest are generated by RDKit. 
    """
    from .utils.validate_and_assign_ids import validate_and_assign

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
    from openff.benchmark.utils.generate_conformers import generate_conformers

    generate_conformers(input_directory,
                        output_directory)
    
if __name__ == "__main__":
    cli()
