"""
Components for energy minimization / geometry optimization of molecules for FF benchmarking.

These components leverage the QCArchive software stack for execution.

"""

import os
import shutil
import logging
logging.disable(logging.WARNING) 

from qcportal import FractalClient

from .seasons import SEASONS

def _mols_from_paths(input_paths):
    from openforcefield.topology import Molecule
    from openforcefield.utils import toolkits
    
    # make sure we deregister OpenEye, if it is present
    try:
        toolkits.GLOBAL_TOOLKIT_REGISTRY.deregister_toolkit(toolkits.OpenEyeToolkitWrapper)
    except toolkits.ToolkitUnavailableException:
        pass

    mols = []
    for path in input_paths:
        if (os.path.isfile(path) and path.split('.')[-1].lower() == 'sdf'):
            mols.append(Molecule.from_file(path))
        elif os.path.isdir(path):
            for root, dirs, files in os.walk(path):
                for file in files:
                    path = os.path.join(root, file)
                    if (os.path.isfile(path) and path.split('.')[-1].lower() == 'sdf'):
                        mols.append(Molecule.from_file(path))

    return mols


def submit_molecules(server_uri, input_paths, season, dataset_name="Benchmark Optimizations", 
        replace=True):
    """Submit SDF molecules from given directory to the target QCFractal server.

    Parameters
    ----------
    server_uri : str
        Target QCFractal server URI.
    input_paths : iterable of Path-like
        Paths to SDF files or directories; if directories, all files SDF files in are loaded, recursively.
    season : str
        Benchmark season identifier. Indicates the mix of compute specs to utilize.
    dataset_name : str
        Dataset name to use for submission on the QCFractal server.

    """
    from qcsubmit.factories import OptimizationDataset, OptimizationDatasetFactory

    # extract molecules from SDF inputs
    mols = _mols_from_paths(input_paths)

    # create OptimizationDataset with QCSubmit
    ds = OptimizationDataset(dataset_name=dataset_name)
    factory = OptimizationDatasetFactory()

    for mol in mols:
        id = "{org}-{molecule:05}-{conformer:02}".format(
                org=mol.properties['group_name'],
                molecule=mol.properties['molecule_index'],
                conformer=mol.properties['conformer_index'])

        attributes = factory.create_cmiles_metadata(mol)
        ds.add_molecule(index=id, molecule=mol, attributes=attributes)

    ds.qc_specifications = SEASONS[season]

    ds.metadata.long_description_url = "https://localhost.local/null"

    logging.info("Submitting...")
    client = FractalClient(server_uri, verify=False)
    ds.submit(client=client)
    logging.info("Submitted!")


def _mol_from_qcserver(record):
    from openforcefield.topology import Molecule
    from openforcefield.utils import toolkits
    
    # make sure we deregister OpenEye, if it is present
    try:
        toolkits.GLOBAL_TOOLKIT_REGISTRY.deregister_toolkit(toolkits.OpenEyeToolkitWrapper)
    except toolkits.ToolkitUnavailableException:
        pass

    return Molecule.from_qcschema(record)


def export_molecule_data(server_uri, output_directory, dataset_name="Benchmark Optimizations",
                         delete_existing=False):
    """Export all molecule data from target QCFractal server to the given directory.

    Parameters
    ----------
    server_uri : str
        Target QCFractal server URI.
    output_directory : str
        Directory path to deposit exported data.
    dataset_name : str
        Dataset name to extract from the QCFractal server.
    delete_existing : bool (False)
        If True, delete existing directory if present.

    """
    from openforcefield.topology.molecule import unit
    import numpy as np

    # get dataset
    client = FractalClient(server_uri, verify=False)
    optds = client.get_collection("OptimizationDataset", dataset_name)
    optds.status()

    try:
        os.makedirs(output_directory)
    except OSError:
        if delete_existing:
            shutil.rmtree(output_directory)
        else:
            raise Exception(f'Output directory {output_directory} already exists. '
                             'Specify `delete_existing=True` to remove.')

    # for each compute spec, create a folder in the output directory
    # deposit SDF giving final molecule, energy
    specs = optds.list_specifications().index.tolist()
    for spec in specs:
        os.makedirs(os.path.join(output_directory, spec))
        optentspec = optds.get_specification(spec)

        records = optds.data.dict()['records']

        for id, opt in optds.df[spec].iteritems():
            mol = _mol_from_qcserver(records[id.lower()])

            # set conformer as final, optimized geometry
            final_qcmol = opt.get_final_molecule()
            geometry = unit.Quantity(
                    np.array(final_qcmol.geometry, np.float), unit.bohr
                )
            mol._add_conformer(geometry.in_units_of(unit.angstrom))

            # fix to ensure output fidelity of ids; losing 02 padding on conformer
            org, molecule, conformer = id.split('-')
            id = "{org}-{molecule:05}-{conformer:02}".format(org=org,
                                                             molecule=int(molecule),
                                                             conformer=int(conformer))

            # set molecule metadata
            mol.name = id

            (mol.properties['group_name'],
             mol.properties['molecule_index'],
             mol.properties['conformer_index']) = id.split('-')

            # SDF key-value pairs should be used for method, basis, program, provenance, `openff-benchmark` version
            mol.properties['method'] = optentspec.qc_spec.method
            mol.properties['basis'] = optentspec.qc_spec.basis
            mol.properties['program'] = optentspec.qc_spec.program

            # SDF key-value pairs also used for final energies
            mol.properties['initial_energy'] = opt.energies[0]
            mol.properties['final_energy'] = opt.energies[-1]

            # subfolders for each compute spec, files named according to molecule ids
            outfile = "{}.sdf".format(
                    os.path.join(output_directory, spec, id))
            mol.to_file(outfile, file_format='sdf')


def get_optimization_status(server_uri, dataset_name="Benchmark Optimizations", client=None):
    """Get status of optimization for each molecule ID

    """
    if client is None:
        client = FractalClient(server_uri, verify=False)

    optds = client.get_collection("OptimizationDataset", dataset_name)
    optds.status()
    return optds.df


def errorcycle_optimizations(server_uri):
    """Error cycle optimizations that have failed.

    """


def execute_optimization_from_server(molecule_id, push_complete=False):
    """Execute optimization from the given molecule locally on this host.

    """
    pass


def execute_optimization_from_molecules(input_paths, output_directory, season,
        ncores=1, delete_existing=False):
    """Execute optimization from the given SDF molecules locally on this host.

    """
    from time import sleep
    from qcfractal import FractalSnowflakeHandler

    # fail early if output_directory already exists and we aren't deleting it
    if os.path.isdir(output_directory):
        if delete_existing:
            shutil.rmtree(output_directory)
        else:
            raise Exception(f'Output directory {output_directory} already exists. '
                             'Specify `delete_existing=True` to remove.')

    dataset_name = "Benchmark Scratch"

    # start up Snowflake
    server = FractalSnowflakeHandler(ncores=ncores)
    client = server.client()
    server_uri = server.get_address()

    # submit molecules
    submit_molecules(server_uri, input_paths, season, dataset_name=dataset_name)

    while True:
        df = get_optimization_status(server_uri, dataset_name, client=client)
        complete = df.applymap(lambda x: x.status.value == 'COMPLETE').sum().sum()
        total = df.size
        print(f"{complete}/{total} COMPLETE", end='\r')
        if complete:
            break
        sleep(5)

    export_molecule_data(server_uri, output_directory, dataset_name=dataset_name)
