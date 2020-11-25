"""
Components for energy minimization / geometry optimization of molecules for FF benchmarking.

These components leverage the QCArchive software stack for execution.

In general, we abstract away much of the server's operational details, including datasets.
We store multiple datasets as needed

"""

import os
import logging

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
    """Submit SDF molecules from given directory.

    Parameters
    ----------
    server_uri : str
        Target QCFractal server URI.
    input_paths : iterable of Path-like
        Paths to SDF files or directories; if directories, all files SDF files in are loaded, recursively.
    season : integer
        Benchmark season number. Indicates the mix of compute specs to utilize.

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


def export_molecule_data(server_uri, destination_path):
    """Export molecule data from target QCFractal instance.

    Parameters
    ----------


    """
    # export all molecule/optimization data from all datasets 


    # SDF key-value pairs should be used for method, basis, program, provenance, `openff-benchmark` version

    # SDF key-value pairs also used for final energies

    # subfolders for each compute spec, files named according to molecule ids


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
