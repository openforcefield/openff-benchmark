"""
Components for energy minimization / geometry optimization of molecules for FF benchmarking.

These components leverage the QCArchive software stack for execution.

"""

import os
import shutil
import logging
import subprocess
from collections import defaultdict

# DO NOT USE; MESSES UP LOGGING 
# consider overriding the loglevel of any particularly noisy loggers individual
#logging.disable(logging.WARNING) 

from ..utils.io import mols_from_paths


def ffbuilder(input_paths, schrodinger_path, host_settings, opls_dir, recursive, output_path):
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
    os.makedirs(f'{output_path}/1-ffbuilder', exist_ok=True)
    
    with open(f'{output_path}/1-ffbuilder/ffb_input.sdf', 'w') as file:
        for mol in mols:
            if mol.name.endswith('00'):
                mol.to_file(file, 'SDF')

    command = [os.path.join(schrodinger_path, 'ffbuilder'),
               '-JOBNAME', 'ffb_openff_benchmark',
               os.path.join(output_path, '1-ffbuilder', 'ffb_input.sdf'),
               '-fit_advanced',
               '-OPLSDIR', opls_dir,
               '-HOST', host_settings,
               '-TMPLAUNCHDIR'
    ]
    subprocess.run(command)
