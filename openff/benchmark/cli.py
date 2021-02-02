"""
Top-level `openff-benchmark` cli entrypoint.

"""

import click
import logging
import csv
import sys

# Deregister OpenEye for any ToolkitRegistry calls that happen in this file
logger = logging.getLogger('openforcefield.utils.toolkits')
prev_log_level = logger.getEffectiveLevel()
logger.setLevel(logging.ERROR)

from openforcefield.utils.toolkits import GLOBAL_TOOLKIT_REGISTRY, OpenEyeToolkitWrapper

logger.setLevel(prev_log_level)

oetk_loaded = False
for tkw in GLOBAL_TOOLKIT_REGISTRY.registered_toolkits:
    if isinstance(tkw, OpenEyeToolkitWrapper):
        oetk_loaded = True
if oetk_loaded:
    GLOBAL_TOOLKIT_REGISTRY.deregister_toolkit(OpenEyeToolkitWrapper)


logging.basicConfig(filename=f'{sys.argv[2]}.log',
                    level=logging.DEBUG
                    )


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

#@optimize.command()
#def submit_compute(fractal_uri, dataset_name, spec_name, method, basis, program):
#    pass

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
@click.option('-f', '--frequency', default=10, help="Time in seconds between updates")
@click.option('-d', '--dataset-name', default=None)
@click.option('-c', '--compute-specs', default=None, help="Comma-separated compute spec names to limit status display to")
def progress(fractal_uri, frequency, dataset_name, compute_specs):
    from time import sleep
    from tqdm import trange

    from .geometry_optimizations.compute import OptimizationExecutor

    if compute_specs is not None:
        compute_specs = compute_specs.split(',')

    optexec = OptimizationExecutor()
    if dataset_name is None:
        datasets = optexec.list_optimization_datasets(fractal_uri=fractal_uri)
    else:
        datasets = [dataset_name]

    dfs = []
    for dataset in datasets:
        dfs.append(optexec.get_optimization_status(fractal_uri=fractal_uri,
                                                   dataset_name=dataset,
                                                   compute_specs=compute_specs))

    complete = sum([len([opt for opt in df.values.flatten() if opt.status == 'COMPLETE']) for df in dfs])
    progbar = trange(sum([df.size for df in dfs]), initial=complete)

    while True:
        dfs = []
        for dataset in datasets:
            dfs.append(optexec.get_optimization_status(
                fractal_uri=fractal_uri, dataset_name=dataset, compute_specs=compute_specs))

        complete_i = sum([len([opt for opt in df.values.flatten() if opt.status == 'COMPLETE']) for df in dfs])
        progbar.update(complete_i - complete)
        progbar.refresh()

        complete = complete_i

        sleep(frequency)


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
    print("\n".join(datasets))


@optimize.command()
@click.option('-u', '--fractal-uri', default="localhost:7777")
@click.argument('dataset-name', nargs=-1)
def delete_datasets(fractal_uri, dataset_name):
    from .geometry_optimizations.compute import OptimizationExecutor

    optexec = OptimizationExecutor()
    datasets = optexec.delete_optimization_datasets(fractal_uri, dataset_name)


@optimize.command()
@click.option('-u', '--fractal-uri', default="localhost:7777")
@click.option('-d', '--dataset-name', required=True)
@click.option('-c', '--compute-specs', default=None, help="Comma-separated compute spec names to limit status display to")
@click.argument('molids', nargs=-1)
def extract(molids, fractal_uri, dataset_name, compute_specs):
    import json
    from .geometry_optimizations.compute import OptimizationExecutor

    optexec = OptimizationExecutor()

    if compute_specs is not None:
        compute_specs = compute_specs.split(',')

    results = optexec.get_optimization_from_server(fractal_uri,
                                                     dataset_name,
                                                     compute_specs=compute_specs,
                                                     molids=molids)

    # export and read back into JSON for final output
    results_processed = [json.loads(res.json()) for res in results]
    print(json.dumps(results_processed))


@optimize.command()
@click.option('-u', '--fractal-uri', default="localhost:7777")
@click.option('-d', '--dataset-name', required=True)
@click.option('-t', '--nthreads', default=2)
@click.option('-m', '--memory', default=2)
@click.option('-c', '--compute-specs', default=None, help="Comma-separated compute spec names to limit status display to")
@click.option('-o', '--output-directory')
@click.option('--scf-maxiter', type=int, default=200)
@click.option('--geometric-maxiter', type=int, default=300)
@click.option('--geometric-coordsys', type=click.Choice(['dlc', 'tric']), default='dlc')
@click.option('--geometric-qccnv', is_flag=True)
@click.argument('molids', nargs=-1)
def execute_from_server(molids, fractal_uri, dataset_name, compute_specs, nthreads, memory, output_directory,
                        scf_maxiter, geometric_maxiter, geometric_coordsys, geometric_qccnv):
    import json
    from .geometry_optimizations.compute import OptimizationExecutor

    optexec = OptimizationExecutor()

    if compute_specs is not None:
        compute_specs = compute_specs.split(',')

    results = optexec.execute_optimization_from_server(fractal_uri,
                                                       dataset_name,
                                                       output_directory=output_directory,
                                                       ncores=nthreads,
                                                       memory=memory,
                                                       compute_specs=compute_specs,
                                                       molids=molids,
                                                       scf_maxiter=scf_maxiter,
                                                       geometric_maxiter=geometric_maxiter,
                                                       geometric_coordsys=geometric_coordsys,
                                                       geometric_qccnv=geometric_qccnv)

    # export and read back into JSON for final output
    results_processed = [json.loads(res.json()) for res in results]
    print(json.dumps(results_processed))


@optimize.command()
@click.option('-s', '--season', required=True)
@click.option('-t', '--nthreads', default=2)
@click.option('-m', '--memory', default=2)
@click.option('--delete-existing', is_flag=True,
              help="Delete existing output directory if it exists")
@click.option('--recursive', is_flag=True, help="Recursively traverse directories for SDF files to include")
@click.option('-o', '--output-directory', required=True, help="Output directory to use for results")
@click.option('--stdout', is_flag=True, help="If set, results also output to STDOUT")
@click.option('--scf-maxiter', type=int, default=200)
@click.option('--geometric-maxiter', type=int, default=300)
@click.option('--geometric-coordsys', type=click.Choice(['dlc', 'tric']), default='dlc')
@click.option('--geometric-qccnv', is_flag=True)
@click.argument('input-paths', nargs=-1)
def execute(input_paths, output_directory, stdout, season, nthreads, memory, delete_existing, recursive,
            scf_maxiter, geometric_maxiter, geometric_coordsys, geometric_qccnv):
    import os
    import json
    import warnings

    from .geometry_optimizations.compute import OptimizationExecutor

    if (not delete_existing) and os.path.exists(output_directory):
        warnings.warn(f"Output directory {output_directory} exists and will be used for export; existing data files will not be replaced.")

    optexec = OptimizationExecutor()

    results = optexec.execute_optimization_from_molecules(input_paths,
                                                          output_directory,
                                                          season,
                                                          ncores=nthreads,
                                                          memory=memory,
                                                          delete_existing=delete_existing,
                                                          recursive=recursive,
                                                          scf_maxiter=scf_maxiter,
                                                          geometric_maxiter=geometric_maxiter,
                                                          geometric_coordsys=geometric_coordsys,
                                                          geometric_qccnv=geometric_qccnv)

    # export and read back into JSON for final output
    if stdout:
        results_processed = [json.loads(res.json()) for res in results]
        print(json.dumps(results_processed))



@cli.group()
def report():
    """Analyze the results and create plots

    """
    pass

@report.command()
@click.option('--input-path', default='./', multiple=True, required=True)
@click.option('--ref-method', default='b3lyp-d3bj', required=True)
@click.option('--output-directory', default='5-compare_forcefields', required=True)
def compare_forcefields(input_path, ref_method, output_directory):
    from .analysis import analysis
    analysis.main(input_path, ref_method, output_directory)

@report.command()
@click.option('--input-path', default='./', multiple=True, required=True)
@click.option('--ref-method', default='b3lyp-d3bj', required=True)
@click.option('--output-directory', default='5-match_minima', required=True)
def match_minima(input_path, ref_method, output_directory):
    from .analysis import analysis
    analysis.match_minima(input_path, ref_method, output_directory)

@report.command()
@click.option('--input-path', default='5-compare_forcefields', multiple=True, required=True)
@click.option('--ref-method', default='default', required=True)
@click.option('--output-directory', default='5-plots-compare-forcefields', required=True)
def plots(input_path, ref_method, output_directory):
    from .analysis import draw
    draw.plot_compare_ffs(input_path, ref_method, output_directory)



@cli.group()
def preprocess():
    """Prepare input molecules for compute.
    """
    pass


@preprocess.command()
@click.option('--delete-existing', is_flag=True)
@click.option('--add',
              is_flag=True,
              help='Appends new molecules to the dataset. The new molecules MUST NOT be new conformers of previously-existing molecules.')
@click.option('-o', '--output-directory',
              default='1-validate_and_assign', 
              help='Directory to put output files. If this directory does not exist, one will be created.')
@click.option('-g', '--group-name',
              required=True,
              help='Group name for assigning IDs to the molecules.')
@click.argument('input-3d-molecules',
                nargs=-1)
def validate(input_3d_molecules, output_directory, group_name, delete_existing, add):
    """
    Validate and assign identifiers to molecules.

    This command preprocesses, validates, and creates a naming system for molecules that will be submitted to the benchmarking workflow. 

    This command takes 3D molecules in SDF format and a three-character "group name" as input. 
    It produces a directory containing:

      * "GGG-MMMMM-CC.sdf": Files containing one validated molecule each, where G is the group name, M is the molecule ID, and C is the conformer index

      * "error_mols": A subdirectory where molecule that do not pass validation are stored

      * "log.txt": Debugging info

      * "name_assignments.csv": A mapping from input file/molecule names to the IDs assigned for this workflow

    At least one 3D geometry of each molecule must be provided.
    The "group name" becomes the first three characters of the validated file names. 
    Multiple conformers of the same molecule will be automatically detected and grouped under a single molecule ID.
    The definition of "identical molecule" is whether RDKit assigns them the same canonical, isomeric, explicit hydrogen SMILES. 
    When molecules are grouped, their RMSD (accounting for symmetry automorphs) is tested. 
    If two inputs are within 1.0 A by RMSD, the second is considered an error and sent to the error_mols directory.
    SD data pairs and molecule names are stripped from the input molecules.
    The order of atoms may change during this step.
    
    
    This command also attempts to detect technical issues that could prevent the files from working with the OpenFF toolkit.
    Files that will cause problems in subsequent steps are routed to the "error_mols" subdirectory of the output directory.
    Where possible, these cases write both an SDF file of the molecule (with key-value pairs indicating the file the structure came from),
    and a correspondingly-named txt file containing more details about the error.
    """
    from openforcefield.topology import Molecule
    from openff.benchmark.utils.validate_and_assign_ids import validate_and_assign
    from openff.benchmark.utils.utils import prepare_folders
    import glob
    import os
    from tqdm import tqdm
    import copy
    import csv

    # If we're not running with the `--add` flag, these will remain as empty lists
    existing_output_mols = []
    existing_name_assignments = []

    # Prepare required directories, ensuring that input flags (`delete_existing` and `add`) are sane
    error_dir = prepare_folders(output_directory=output_directory, delete_existing=delete_existing, add=add)

    input_mols = []
    error_mols = []
    # Load all input mols, annotating SD data with original file+mol index
    print('Loading input files')
    logging.info('Loading input files')
    for molecule_3d_file in tqdm(input_3d_molecules):
        assert os.path.exists(molecule_3d_file), f"File {molecule_3d_file} does not exist"
        # Squelch reading warnings
        # We'll recover the full text of the warning/error during the file round trip tests in validate_and_assign_ids.py.
        try:
            toolkit_logger = logging.getLogger('openforcefield.utils.toolkits')
            prev_log_level = toolkit_logger.getEffectiveLevel()
            toolkit_logger.setLevel(logging.ERROR)
            loaded_mols = Molecule.from_file(molecule_3d_file,
                                             file_format='sdf',
                                             allow_undefined_stereo=True)
            toolkit_logger.setLevel(prev_log_level)
        except Exception as e:
            error_mols.append((molecule_3d_file, None, e))
            continue

        if not isinstance(loaded_mols, list):
            loaded_mols = [loaded_mols]
        for mol_index, mol in enumerate(loaded_mols):
            # Append the most recent file info to the provenance properties
            # These properties will be in the same order as the molecule's conformers.
            mol.properties['original_file'] = molecule_3d_file
            mol.properties['original_file_index'] = mol_index
            mol.properties['original_name'] =mol.name
            input_mols.append(copy.deepcopy(mol))
    logging.info('input_mols')
    logging.info(input_mols)
    # Load all pre-existing output mols
    if add:
        print('Running in `--add` mode. Loading existing output molecules.')
        logging.info('Running in `--add` mode. Loading existing output molecules.')
        existing_output_mol_files = glob.glob(os.path.join(output_directory, '*.sdf'))
        existing_output_mols = []
        for existing_file in tqdm(existing_output_mol_files):
            existing_output_mols.append(Molecule.from_file(existing_file))
        # Strip the "-CC.sdf" from the output names ("GGG-MMMMM-CC.sdf")
        existing_name_assignments = []
        with open(os.path.join(output_directory, 'name_assignments.csv')) as of:
            csv_reader = csv.reader(of)
            for row in csv_reader:
                existing_name_assignments.append(row)
            # Remove header line
            existing_name_assignments = existing_name_assignments[1:]

    output = validate_and_assign(input_mols,
                                 group_name,
                                 add,
                                 existing_output_mols,
                                 existing_name_assignments)

    success_mols, new_error_mols, name_assignments = output
    error_mols += new_error_mols
    logging.info('output')
    logging.info(output)
    logging.info('error_mols')
    logging.info(error_mols)
    # Write successfully processed mols
    print('Writing successful mols')
    logging.info('Writing successful mols')
    for success_mol in tqdm(success_mols):
        success_mol.to_file(os.path.join(output_directory, success_mol.name + ".sdf"), "sdf")


    # Write error mols
    print('Writing error mols (if any)')
    logging.info('Writing error mols (if any)')

    # In case error molecules already exist, make sure we don't overwrite them
    existing_error_mols = glob.glob(os.path.join(error_dir, '*sdf'))
    existing_error_mols = [os.path.basename(i) for i in existing_error_mols]
    if len(existing_error_mols) == 0:
        start_error_index = 0
    else:
        existing_error_mol_indices = [int(i.replace('error_mol_','').replace('.sdf', '')) for i in existing_error_mols]
        start_error_index = max(existing_error_mol_indices) + 1
    for idx, (filename, error_mol, exception) in enumerate(tqdm(error_mols)):
        unique_idx = idx + start_error_index
        output_mol_file = os.path.join(error_dir, f'error_mol_{unique_idx}.sdf')
        try:
            error_mol.to_file(output_mol_file, file_format='sdf')
        except Exception as e:
            exception = str(exception)
            exception += "\n Then failed when trying to write mol to error directory with "
            exception += str(e)
        output_summary_file = os.path.join(error_dir, f'error_mol_{unique_idx}.txt')
        with open(output_summary_file, 'w') as of:
            of.write(f'source: {filename}\n')
            of.write(f'error text: {exception}\n')

    # Write name assignments
    import csv
    print('Writing name assignments')
    logging.info('Writing name assignments')
    with open(os.path.join(output_directory, 'name_assignments.csv'), 'w', newline='') as of:
        of.write('orig_name,orig_file,orig_file_index,out_file_name\n')
        csvw = csv.writer(of)
        #np.savetxt(out_file_name, np.array(name_assignments))
        for name_assignment in name_assignments:
            csvw.writerow(name_assignment)



@preprocess.command()
@click.option('--delete-existing', is_flag=True)
@click.option('--add',
              is_flag=True,
              help='Appends new molecules to the dataset. The new molecules MUST be new unique molecules, NOT new conformers of previously-existing molecules.')
@click.option('-o', '--output-directory',
              default='2-generate_conformers', 
              help='Directory for output files. If this directory does not exist, one will be created.')
@click.argument('input-directory')
def generate_conformers(input_directory, output_directory, add, delete_existing):
    """Generate additional conformers for validated molecules.

    INPUT-DIRECTORY should contain validated molecules;
    this is the output directory from the validate step.

    Up to ten distinct conformers of each molecule will be output from this step.
    User-provided conformers from the previous step are always in this output, and have the same molecule and conformer IDs. 
    The conformers generated by this step have at least a 1.5 A heavy-atom RMSD from each other and the already-existing conformers.
    If a 1.5 A heavy atom cutoff produces more than 10 conformers, then the cutoff is gradually incremented until 10 or less conformers remain.
    The conformer pruning uses a "greedy" mechanism, incrementing upwards through conformer ID and pruning all higher-numbered conformers within the cutoff. 
    

    """
    from openforcefield.topology import Molecule
    from openff.benchmark.utils.generate_conformers import generate_conformers
    from openff.benchmark.utils.utils import prepare_folders
    import glob
    import os
    import shutil
    import re


    # Prepare required directories, ensuring that input flags (`delete_existing` and `add`) are sane
    error_dir = prepare_folders(output_directory=output_directory, delete_existing=delete_existing, add=add)

    input_files = glob.glob(os.path.join(input_directory,'*.sdf'))

    if add:
        # now we need to load output molecules and get the difference between them
        output_files = set([os.path.split(path)[-1] for path in  glob.glob(os.path.join(output_directory, "*00.sdf"))])
        # make a new input file list
        in_file_basenames = set([os.path.split(path)[-1] for path in input_files])
        # now get the difference between them
        new_file_basenames = set(in_file_basenames) - output_files
        # if no new files exit here
        if not new_file_basenames:
            print(f"No new files found in {input_directory}, the coverage report was not changed.")
            return
        # now remake the file paths
        molecule_files = [os.path.join(input_directory, filename) for filename in new_file_basenames]

    else:
        molecule_files = input_files
    group2idx2mols2confs = {}

    for input_file in molecule_files:
        # try:
        # We have to allow undefined stereo here because of some silliness with
        # RDKitToolkitWrapper percieveing stereo around terminal ethene groups
        # Note: Allowing undefined stereo shouldn't be necessary any more, ever since the 0.8.2 release
        # See https://github.com/openforcefield/openff-toolkit/pull/786
        mol = Molecule.from_file(input_file, allow_undefined_stereo=True)
        # except:
        #    raise Exception(input_3d_file)
        match = re.findall('([a-zA-Z0-9]{3})-([0-9]{5})-([0-9]{2}).sdf',
                           input_file)
        assert len(match) == 1
        group_id = match[0][0]
        mol_idx = match[0][1]
        conf_id = match[0][2]
        if group_id not in group2idx2mols2confs:
            group2idx2mols2confs[group_id] = {}
        if mol_idx not in group2idx2mols2confs[group_id]:
            group2idx2mols2confs[group_id][mol_idx] = {}
        assert conf_id not in group2idx2mols2confs[group_id][mol_idx]

        group2idx2mols2confs[group_id][mol_idx][conf_id] = mol
        # raise Exception((group_id, mol_id))

    # Call the underlying method to perform conformer generation
    outputs = generate_conformers(group2idx2mols2confs)
    success_mols, error_mols = outputs

    # Write successfully processed mols
    for success_mol in success_mols:
        success_mol.to_file(os.path.join(output_directory, success_mol.name + ".sdf"), "sdf")


    # Write error mols
    error_dir = os.path.join(output_directory, 'error_mols')
    for error_mol, exception in error_mols:
        output_mol_file = os.path.join(error_dir, f'{error_mol.name}.sdf')
        try:
            error_mol.to_file(output_mol_file, file_format='sdf')
        except Exception as e:
            exception = str(exception)
            exception += "\n Then failed when trying to write mol to error directory with "
            exception += str(e)
        output_summary_file = os.path.join(error_dir, f'{error_mol.name}.txt')
        with open(output_summary_file, 'w') as of:
            of.write(f'source: {error_mol.name}\n')
            of.write(f'error text: {exception}\n')


@preprocess.command()
@click.argument("input_directory")
@click.option('--add',
              is_flag=True,
              help='Appends new molecules to the dataset. The new molecules MUST be new unique molecules, NOT new conformers of previously-existing molecules.')
@click.option("-ff", "--forcefield_name", default="openff_unconstrained-1.3.0.offxml")
@click.option("-o", "--output_directory",
              default="3-coverage_report",
              help="The directory for output files.")
@click.option("-p", "--processors",
              default=None,
              type=click.INT)
@click.option('--delete-existing', is_flag=True)
def coverage_report(input_directory, forcefield_name, output_directory, processors, delete_existing, add):
    """
    Generate a coverage report for the set of validated input molecules.
    """
    from openff.benchmark.utils.coverage_report import generate_coverage_report, _update_coverage
    from openforcefield.topology import Molecule
    from openff.benchmark.utils.utils import prepare_folders
    import json
    import glob
    import os
    import shutil

    error_dir = prepare_folders(output_directory=output_directory, delete_existing=delete_existing, add=add)
    # Search for the 00th conformer so we dont double-count any moleucles
    input_files = glob.glob(os.path.join(input_directory, "*00.sdf"))

    if add:
        # now we need to load output molecules and get the difference between them
        output_files = set([os.path.split(path)[-1] for path in  glob.glob(os.path.join(output_directory, "*00.sdf"))])
        # make a new input file list
        in_files = set([os.path.split(path)[-1] for path in input_files])
        # now get the difference between them
        new_files = set(in_files) - output_files
        # if no new files exit here
        if not new_files:
            print(f"No new files found in {input_directory}, the coverage report was not changed.")
            return
        # now remake the file paths
        molecule_files = [os.path.join(input_directory, filename) for filename in new_files]

    else:
        molecule_files = input_files

    # now load each molecule they should already be unique
    molecules = [Molecule.from_file(mol_file, file_format="sdf", allow_undefined_stereo=True) for mol_file in molecule_files]

    report, success_mols, error_mols = generate_coverage_report(input_molecules=molecules,
                                      forcefield_name=forcefield_name,
                                      processors=processors)

    # Copy successfully processed mols and all conformers to the new folder
    for success_mol in success_mols:
        common_id = f"{success_mol.properties['group_name']}-{success_mol.properties['molecule_index']}"
        # get all conformer files
        conformer_files = glob.glob(os.path.join(input_directory, f"{common_id}*"))
        for file in conformer_files:
            shutil.copy(file, output_directory)
    # Copy all conformers of an error mol to the error dir
    for error_mol, e in error_mols:
        with open(os.path.join(error_dir, error_mol.name + ".txt"), 'w') as of:
            of.write(f'source: {error_mol.name}\n')
            of.write(f'error text: {e}\n')
        # now get all conformers and move them
        common_id = f"{error_mol.properties['group_name']}-{error_mol.properties['molecule_index']}"
        conformer_files = glob.glob(os.path.join(input_directory, f"{common_id}*"))
        for file in conformer_files:
            shutil.copy(file, error_dir)

    # Write coverage report
    if add:
        # load the old coverage report
        with open(os.path.join(output_directory, "coverage_report.json")) as old_data:
            old_report = json.load(old_data)

        # now extend the old report with new values
        total_molecules = old_report.pop("total_unique_molecules") + report.pop("total_unique_molecules")
        passed_molecules = old_report.pop("passed_unique_molecules") + report.pop("passed_unique_molecules")
        old_forcefield = old_report.pop("forcefield_name")
        new_forcefield = report.pop("forcefield_name")
        assert old_forcefield == new_forcefield
        # now add
        _update_coverage(old_report, report)
        # now put the totals back in
        old_report["total_unique_molecules"] = total_molecules
        old_report["passed_unique_molecules"] = passed_molecules
        old_report["forcefield_name"] = old_forcefield
        data = json.dumps(old_report, indent=2)

    else:
        data = json.dumps(report, indent=2)

    with open(os.path.join(output_directory, "coverage_report.json"), "w") as reporter:
        reporter.write(data)
        # TODO do we want the list of errors in the coverage report as well?
    

if __name__ == "__main__":
    cli()
