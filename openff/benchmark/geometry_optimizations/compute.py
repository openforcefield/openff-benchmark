"""
Components for energy minimization / geometry optimization of molecules for FF benchmarking.

These components leverage the QCArchive software stack for execution.

"""

import os
import logging
logging.disable(logging.WARNING) 

from openforcefield.topology import Molecule
from openforcefield.utils import toolkits

# make sure we deregister OpenEye, if it is present
try:
    toolkits.GLOBAL_TOOLKIT_REGISTRY.deregister_toolkit(toolkits.OpenEyeToolkitWrapper)
except toolkits.ToolkitUnavailableException:
    pass

from qcsubmit.factories import OptimizationDataset, OptimizationDatasetFactory
from qcportal import FractalClient

from .seasons import SEASONS


def submit_molecules(server_uri, input_paths, season, dataset_name="Benchmark Optimizations"):
    """Submit SDF molecules from given directory to the target QCFractal server.

    Parameters
    ----------
    server_uri : str
        Target QCFractal server URI.
    input_paths : iterable of Path-like
        Paths to SDF files or directories; if directories, all files SDF files in are loaded, recursively.
    season : integer
        Benchmark season number. Indicates the mix of compute specs to utilize.
    dataset_name : str
        Dataset name to use for submission on the QCFractal server.

    """
    # extract molecules from SDF inputs
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

    # create OptimizationDataset with QCSubmit
    ds = OptimizationDataset(dataset_name=dataset_name)
    factory = OptimizationDatasetFactory()

    for mol in mols:
        attributes = factory.create_cmiles_metadata(mol)
        ds.add_molecule(index=mol.properties['id'], molecule=mol, attributes=attributes)

    ds.qc_specifications = SEASONS[season]

    ds.metadata.long_description_url = "https://localhost.local/null"

    logging.info("Submitting...")
    client = FractalClient(server_uri, verify=False)
    ds.submit(client=client)
    logging.info("Submitted!")


def export_molecule_data(server_uri, output_directory, dataset_name="Benchmark Optimizations"):
    """Export all molecule data from target QCFractal server to the given directory.

    Parameters
    ----------
    server_uri : str
        Target QCFractal server URI.
    destination_path : str
        Directory path to deposit exported data.
    dataset_name : str
        Dataset name to extract from the QCFractal server.

    """
    # get dataset
    client = FractalClient(server_uri, verify=False)
    optds = client.get_collection("OptimizationDataset", dataset_name)
    optds.status()

    try:
        os.makedirs(output_directory)
    except OSError:
        raise Exception(f'Output directory {output_directory} already exists. '
                         'The user must delete this manually.')

    # for each compute spec, create a folder in the output directory
    # deposit SDF giving final molecule, energy
    specs = optds.list_specifications().index.tolist()
    for spec in specs:
        os.makedirs(os.path.join(output_directory, spec))
        optentspec = optds.get_specification(spec)

        for id, opt in optds.df[spec].iteritems():
            mol = Molecule.from_qcschema(
                    optds.data.dict()['records']['off-00000-0'], client=client)
            #qcmol_d = opt.get_final_molecule().dict(encoding='json')
            #qcmol_d['attributes'] = qcmol_d['extras']
            #mol = Molecule.from_qcschema(qcmol_d, client=client)
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


def get_optimization_status(server_uri):
    """Get status of optimization for each molecule ID

    """


def errorcycle_optimizations(server_uri):
    """Error cycle optimizations that have failed.

    """


def execute_optimization_from_server(molecule_id, push_complete=False):
    """Execute optimization from the given molecule locally on this host.

    """


def execute_optimization_from_molecule(molecule_path):
    """Execute optimization from the given molecule locally on this host.

    """
