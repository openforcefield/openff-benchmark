"""
Top-level `openff-benchmark` cli entrypoint.

"""

import click

@click.group()
def cli():
    pass


@cli.group()
def optimize():
    """Execute benchmarking geometry optimizations.

    """
    pass

@optimize.command()
@click.option('-u', '--fractal-uri', default="localhost:7777")
@click.option('-d', '--dataset-name', required=True, help="Dataset name to submit molecules under")
@click.option('--recursive', is_flag=True, help="Recursively traverse directories for SDF files to submit")
@click.option('-s', '--season', required=True)
@click.argument('input-path', nargs=-1)
def submit_molecules(fractal_uri, input_path, season, dataset_name, recursive):
    from .geometry_optimizations.compute import OptimizationExecutor

    optexec = OptimizationExecutor()
    optexec.submit_molecules(
            fractal_uri, input_path, season, dataset_name, recursive=recursive)

@optimize.command()
def submit_compute(fractal_uri, dataset_name, spec_name, method, basis, program):
    pass

@optimize.command()
@click.option('-u', '--fractal-uri', default="localhost:7777")
@click.option('-d', '--dataset-name', required=True)
@click.option('--delete-existing', is_flag=True)
@click.option('-o', '--output-directory', default='3-export-compute')
def export(fractal_uri, output_directory, dataset_name, delete_existing):
    import os
    import warnings
    from .geometry_optimizations.compute import OptimizationExecutor

    if (not delete_existing) and os.path.exists(output_directory):
        warnings.warn(f"Output directory {output_directory} exists and will be used for export; existing data files will not be replaced.")

    optexec = OptimizationExecutor()
    optexec.export_molecule_data(
            fractal_uri, output_directory, dataset_name, delete_existing)

@optimize.command()
@click.option('-u', '--fractal-uri', default="localhost:7777")
@click.option('-d', '--dataset-name', required=True)
@click.option('-c', '--compute-specs', default=None, help="Comma-separated compute spec names to limit status display to")
@click.argument('molids', nargs=-1)
def status(molids, fractal_uri, dataset_name, compute_specs):
    import pandas as pd
    from .geometry_optimizations.compute import OptimizationExecutor

    optexec = OptimizationExecutor()

    if compute_specs is not None:
        compute_specs = compute_specs.split(',')

    optdf = optexec.get_optimization_status(fractal_uri, dataset_name, compute_specs=compute_specs, molids=molids)
    pd.options.display.max_rows = len(optdf)
    print(optdf.applymap(lambda x: x.status.value))


@optimize.command()
@click.option('-u', '--fractal-uri', default="localhost:7777")
@click.option('-d', '--dataset-name', required=True)
@click.argument('priority', nargs=1)
def set_priority(fractal_uri, dataset_name, priority):
    from .geometry_optimizations.compute import OptimizationExecutor

    if priority not in ('high', 'normal', 'low'):
        raise ValueError("PRIORITY must be one of 'high', 'normal', or 'low'")

    optexec = OptimizationExecutor()
    optexec.set_optimization_priority(fractal_uri, priority, dataset_name)


@optimize.command()
@click.option('-u', '--fractal-uri', default="localhost:7777")
@click.option('-d', '--dataset-name', required=True)
@click.option('-c', '--compute-specs', default=None, help="Comma-separated compute spec names to limit errorcycling to")
@click.option('-w', '--watch', default=None, help="Run in 'watch' mode with interval in seconds between cycles")
@click.argument('molids', nargs=-1)
def errorcycle(molids, fractal_uri, dataset_name, watch, compute_specs):
    from time import sleep
    from .geometry_optimizations.compute import OptimizationExecutor

    optexec = OptimizationExecutor()

    if compute_specs is not None:
        compute_specs = compute_specs.split(',')

    if watch is not None:
        while True:
            optexec.errorcycle_optimizations(fractal_uri, dataset_name, compute_specs=compute_specs, molids=molids)
            sleep(watch)
    else:
        optexec.errorcycle_optimizations(fractal_uri, dataset_name, compute_specs=compute_specs, molids=molids)


@optimize.command()
@click.option('-u', '--fractal-uri', default="localhost:7777")
@click.option('-d', '--dataset-name', required=True)
@click.option('-c', '--compute-specs', default=None, help="Comma-separated compute spec names to limit status display to")
@click.argument('molids', nargs=-1)
def traceback(molids, fractal_uri, dataset_name, compute_specs):
    import pandas as pd
    from .geometry_optimizations.compute import OptimizationExecutor

    optexec = OptimizationExecutor()

    if compute_specs is not None:
        compute_specs = compute_specs.split(',')

    tracebacks = optexec.get_optimization_tracebacks(fractal_uri,
                                                     dataset_name, 
                                                     compute_specs=compute_specs, 
                                                     molids=molids)
    for index, row in tracebacks.iterrows():
        print('\n' + index)
        for column in row.index:
            print('-' * len(index))
            print('#', column, '\n')
            print(row[column])
        print('==============')


@optimize.command()
@click.option('-u', '--fractal-uri', default="localhost:7777")
def list_datasets(fractal_uri):
    import pandas as pd
    from .geometry_optimizations.compute import OptimizationExecutor

    optexec = OptimizationExecutor()

    datasets = optexec.list_optimization_datasets(fractal_uri)
    print("\n".join(datasets.reset_index()['name'].tolist()))

@optimize.command()
@click.option('-u', '--fractal-uri', default="localhost:7777")
@click.argument('dataset-name', nargs=-1)
def delete_datasets(fractal_uri, dataset_name):
    from .geometry_optimizations.compute import OptimizationExecutor

    optexec = OptimizationExecutor()
    datasets = optexec.delete_optimization_datasets(fractal_uri, dataset_name)

#@optimize.command()
#def execute_from_server():
#    pass

#@optimize.command()
#def debug_from_server():
#    pass

@optimize.command()
@click.option('-s', '--season', required=True)
@click.option('--ncores', default=2)
@click.option('--delete-existing', is_flag=True,
              help="Delete existing output directory if it exists")
@click.option('--recursive', is_flag=True, help="Recursively traverse directories for SDF files to include")
@click.option('-o', '--output-directory', default='3-export-compute')
@click.argument('input-paths', nargs=-1)
def execute(input_paths, output_directory, season, ncores, delete_existing, recursive):
    import os
    import signal
    import warnings

    from .geometry_optimizations.compute import OptimizationExecutor

    if (not delete_existing) and os.path.exists(output_directory):
        warnings.warn(f"Output directory {output_directory} exists and will be used for export; existing data files will not be replaced.")

    optexec = OptimizationExecutor()
    dataset_name='Benchmark Scratch'

    def handle_signal(sig, frame):
        optexec.stop = True

    signal.signal(signal.SIGINT, handle_signal)
    signal.signal(signal.SIGTERM, handle_signal)

    optexec.execute_optimization_from_molecules(
            input_paths, output_directory, season, ncores=ncores,
            delete_existing=delete_existing,
            recursive=recursive)



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
@click.option('--ref_method', default='default', required=True)
@click.option('--output_directory', default='./results', required=True)
def match_minima(input_path, ref_method, output_directory):
    from .analysis import analysis
    analysis.match_minima(input_path, ref_method, output_directory)

@report.command()
@click.option('--input-path', default='./', multiple=True, required=True)
def plots(input_path):
    from .analysis import draw
    draw.plot_compare_ffs(input_path)



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
