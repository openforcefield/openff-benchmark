"""
Components for energy minimization / geometry optimization of molecules for FF benchmarking.

These components leverage the QCArchive software stack for execution.

In general, we abstract away much of the server's operational details, including datasets.
We store multiple datasets as needed

"""

def submit_molecules(server_uri, input_path, compute_specs):
    """Submit SDF molecules from given directory.

    """
    # assume files in a single directory

    # if a filepath is given, submit only that molecule

    # each time this is called, a new dataset is created on the server
    # this is an implementation detail


def export_molecule_data(server_uri, destination_path):
    """Export molecule data from target QCFractal instance.

    """
    # export all molecule/optimization data from all datasets 


    # SDF key-value pairs should be used for method, basis, program, provenance, `openff-benchmark` version

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
