"""
Components for energy minimization / geometry optimization of molecules for FF benchmarking.

These components leverage the QCArchive software stack for execution.

"""

import os
import shutil
import logging
import subprocess
import warnings
from collections import defaultdict

# DO NOT USE; MESSES UP LOGGING
# consider overriding the loglevel of any particularly noisy loggers individual
# logging.disable(logging.WARNING)

from ..utils.io import mols_from_paths


def ffbuilder(
    input_paths, schrodinger_path, host_settings, opls_dir, recursive, output_path
):
    """Submit SDF molecules from given directory to a Schrodinger ffbuilder OPLS3e parameterization

    Parameters
    ----------
    input_paths : str or list of str
        Input SDF molecules or paths containing SDFs
    schrodinger_path : str
        Path to Schrodinger executables
    host_settings : str
        Settings for Schrodinger job of format '<host>:<max number of jobs>'
    opls_dir : str
        Path to custom OPLS3e parameter set
    recursive : bool
        If True, recursively load SDFs from any directories given in `input_paths`.
    output_path: str
        Path where output sdf files should be written to
    """
    # extract molecules from SDF inputs
    mols = mols_from_paths(input_paths, recursive=recursive)

    os.makedirs(output_path, exist_ok=True)
    os.makedirs(f"{output_path}", exist_ok=True)

    if os.path.exists(os.path.join(output_path, "ffb_openff_benchmark.zip")):
        warnings.warn(
            f"Schrodinger ffbuilder has been run already. "
            f"Merging parameters and moving files to backup."
        )
        ffb_merge(schrodinger_path, opls_dir, output_path)
        ffb_backup(output_path)

    mol_ids = set()
    with open(f"{output_path}/ffb_input.sdf", "w") as file:
        for mol in mols:
            mol_id = mol.properties["molecule_index"]
            if mol_id in mol_ids:
                pass
            else:
                mol_ids.add(mol_id)
                mol.to_file(file, "SDF")
    if len(mol_ids) == 0:
        raise Exception("No molecules are selected for the ffbuilder job. "
                        "Check the input path.")

    command = [
        os.path.join(schrodinger_path, "ffbuilder"),
        "-JOBNAME",
        "ffb_openff_benchmark",
        "ffb_input.sdf",
        "-OPLSDIR",
        opls_dir,
        "-HOST",
        host_settings,
        "-TMPLAUNCHDIR",
    ]
    subprocess.run(command, cwd=output_path)


def ffb_merge(schrodinger_path, opls_dir, input_path):
    if os.path.exists(os.path.join(input_path, "ffb_openff_benchmark_oplsdir")):
        command = [
            os.path.join(schrodinger_path, "run"),
            "-FROM",
            "ffld",
            "merge_opls_data.py",
            os.path.join(input_path, "ffb_openff_benchmark_oplsdir"),
            "-o",
            opls_dir,
        ]
        subprocess.run(command)


def ffb_backup(output_path):
    os.makedirs(f"{output_path}", exist_ok=True)
    backup_base = "ffb_openff_benchmark.bk"
    backup_dir = backup_base
    num = 1
    while os.path.exists(os.path.join(output_path, backup_dir)):
        backup_dir = f"{backup_base}.{num}"
        num += 1
    backup_path = os.path.join(f"{output_path}", backup_dir)
    os.makedirs(backup_path)
    for obj in [
        "ffb_openff_benchmark_oplsdir",
        "ffb_input.sdf",
        "ffb_openff_benchmark.log",
        "ffb_openff_benchmark.zip",
    ]:
        obj_path = os.path.join(f"{output_path}", obj)
        if os.path.exists(obj_path):
            shutil.move(obj_path, backup_path)
