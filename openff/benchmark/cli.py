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
    """
    Validate and assign identifiers to molecules.

    This command preprocesses, validates, and creates a naming system for molecules that will be submitted to the benchmarking workflow. 

    This command takes 3D molecules in SDF format and a three-character "group name" as input. 
    It produces a directory containing:

      * "GGG-MMMMM-CC.sdf": Files containing one validated molecule each, where G is the group name, M is the molecule ID, and C is the conformer index

      * "error_mols": A directory where molecule that do not pass validation are stored

      * "log.txt": Debugging info

      * "name_assignments.csv": A mapping from input file/molecule names to the IDs assigned for this workflow

    At least one 3D geometry of each molecule must be provided.
    The "group name" becomes the first three characters of the validated file names. 
    Multiple confomers of the same molecule will be automatically detected and grouped under a single molecule ID.
    The definition of "identical molecule" is whether RDKit assigns them the same canonical, isomeric, explicit hydrogen SMILES. 
    When molecules are grouped, their RMSD (accounting for symmetry automorphs) is tested. 
    If two inputs are within 0.1 A by RMSD, the second is considered an error and sent to the error_mols directory.
    SD data pairs and molecule names are stripped from the input molecules.
    The order of atoms may change during this step.
    
    
    This command also attempts to detect technical issues that could prevent the files from working with the OpenFF toolkit. Files that will cause problems in subsequent steps are routed to the "error_mols" subdirectory of the output directory. Where possible, these cases write both an SDF file of the molecule (with key-value paris indicating the file the structure came from), and a correspondingly-named txt file containing more details about the error.
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

    Up to ten distinct conformers of each molecule will be output from this step.
    User-provided conformers from the previous step are always in this output, and have the same molecule and conformer IDs. 
    The conformers generated by this step have at least a 1 A heavy-atom RMSD from each other.
    If a 1 A heavy atom cutoff produces more than 10 conformers, then the cutoff is gradually incremented until 10 or less conformers remain.
    The conformer pruning uses a "greedy" mechanism, incrementing upwards through conformer ID and pruning all higher-numbered conformers within the cutoff. 
    

    """
    from openff.benchmark.utils.generate_conformers import generate_conformers

    generate_conformers(input_directory,
                        output_directory,
                        delete_existing=delete_existing)
    
if __name__ == "__main__":
    cli()
