"""
Top-level `openff-benchmark` cli entrypoint.

"""

import click

@click.group()
def cli():
    pass


@cli.group()
def fractal():
    """Initialize and start a QCFractal server and managers.

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
    from .geometry_optimizations.compute import OptimizationExecutor

    optexec = OptimizationExecutor()
    optexec.submit_molecules(
            server_uri, input_path, season, dataset_name)

@optimize.command()
def submit_compute(server_uri, input_path):
    pass

@optimize.command()
@click.option('--server-uri', default="localhost:7777")
@click.option('--dataset-name', required=True)
@click.option('--delete-existing', is_flag=True)
@click.option('--keep-existing', is_flag=True,
              help="Keep using existing output directory if it exists, but do not rewrite results for IDs already present")
@click.option('-o', '--output-directory', default='3-export-compute')
def export(server_uri, output_directory, dataset_name, delete_existing, keep_existing):
    import os
    from .geometry_optimizations.compute import OptimizationExecutor

    if keep_existing and delete_existing:
        raise ValueError("Cannot use both `--delete-existing` and `--keep-existing`; choose one or neither.")

    if not (delete_existing or keep_existing) and os.path.exists(output_directory):
        raise ValueError(f"Output directory {output_directory} exists; specify `--delete-existing` `--keep-existing` to proceed.")

    optexec = OptimizationExecutor()
    optexec.export_molecule_data(
            server_uri, output_directory, dataset_name, delete_existing, keep_existing)

@optimize.command()
@click.option('--server-uri', default="localhost:7777")
@click.option('--dataset-name', required=True)
def status(server_uri, dataset_name):
    from .geometry_optimizations.compute import OptimizationExecutor

    optexec = OptimizationExecutor()
    optdf = optexec.get_optimization_status(server_uri, dataset_name)
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
@click.option('--delete-existing', is_flag=True,
              help="Delete existing output directory if it exists")
@click.option('--keep-existing', is_flag=True,
              help="Keep using existing output directory if it exists, but do not recalculate results for IDs already present")
@click.option('-o', '--output-directory', default='3-export-compute')
@click.argument('input-path', nargs=-1)
def execute(input_path, output_directory, season, ncores, delete_existing, keep_existing):
    import signal

    from qcfractal import FractalSnowflakeHandler

    from .geometry_optimizations.compute import OptimizationExecutor

    optexec = OptimizationExecutor()

    if keep_existing and delete_existing:
        raise ValueError("Cannot use both `--delete-existing` and `--keep-existing`; choose one or neither.")

    # start up Snowflake
    server = FractalSnowflakeHandler(ncores=ncores)
    dataset_name='Benchmark Scratch'

    def handle_signal(sig, frame):
        print('stopping!')
        optexec.stop = True
        #server_uri = server.get_address()

        ## one final export of data
        #optcompute.export_molecule_data(server_uri, output_directory,
        #        dataset_name=dataset_name, delete_existing=False, keep_existing=True)

        ## stop server
        #server.stop()

    signal.signal(signal.SIGINT, handle_signal)
    signal.signal(signal.SIGTERM, handle_signal)

    optexec.execute_optimization_from_molecules(
            server, input_path, output_directory, season,
            delete_existing=delete_existing, keep_existing=keep_existing)

    print('stopping server')
    server.stop()


@cli.group()
def preprocess():
    """Prepare input molecules for compute.

    """
    pass

@preprocess.command()
@click.option('--delete-existing', is_flag=True)
@click.option('-o', '--output-directory',
              default='1-validate_and_assign', 
              help='Directory to put output files. If this directory does not exist, one will be created.')
@click.option('-g', '--group-name',
              required=True,
              help='Group name for assigning IDs to the molecules.')
@click.argument('input-3d-molecules',
                nargs=-1)
def validate(input_3d_molecules, output_directory, group_name, delete_existing):
    """Validate and assign identifiers to molecules.

    This script preprocesses, validates, and creates a naming system for molecles that will be submitted to the benchmarking workflow. 
    For each unique molecule, up to ten total conformers will be generated by the benchmarking workflow.\n
    This script takes two forms of molecule inputs:\n
        1) "Graph Molecules", the "normal" input type where the geometry of the input molecule will NOT be passed to subsequent steps, and instead all conformers of the molecule will be generated by RDKit.\n
        2) "3D Molecules", where one or more specific geometries are provided, and the rest are generated by RDKit. 
    """
    from .utils.validate_and_assign_ids import validate_and_assign

    validate_and_assign(input_3d_molecules,
                        group_name,
                        output_directory,
                        delete_existing=delete_existing)


@preprocess.command()
@click.option('--delete-existing', is_flag=True)
@click.option('-o', '--output-directory',
              default='2-generate_conformers', 
              help='Directory for output files. If this directory does not exist, one will be created.')
@click.argument('input-directory')
def generate_conformers(input_directory, output_directory, delete_existing):
    """Generate additional conformers for validated molecules.

    INPUT-DIRECTORY should contain validated molecules;
    this is the output directory from the validate step.

    """
    from openff.benchmark.utils.generate_conformers import generate_conformers

    generate_conformers(input_directory,
                        output_directory,
                        delete_existing=delete_existing)
    
if __name__ == "__main__":
    cli()
