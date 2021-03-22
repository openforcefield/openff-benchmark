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

from openforcefield.topology import Molecule
from ..utils.io import mols_from_paths

def sdf_to_mae(sdf_file, mae_file, schrodinger_path, append=True):
    command = [os.path.join(schrodinger_path, 'utilities', 'sdconvert'),
               '-isd', sdf_file,
               '-omae', mae_file,
    ]
    if append:
        command.append('-a')
    subprocess.run(command)

def mae_to_sdf(mae_file, sdf_file, schrodinger_path, append=True):
    command = [os.path.join(schrodinger_path, 'utilities', 'sdconvert'),
               '-imae', mae_file,
               '-osd', sdf_file,
    ]
    if append:
        command.append('-a')
    subprocess.run(command)

def optimization(input_paths, schrodinger_path, host_settings, opls_dir, recursive, output_path):
    """Submit SDF molecules from given directory to a Schrodinger macromodel OPLS3e minimization

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
    os.makedirs(output_path, exist_ok=True)
    mae_file = os.path.join(output_path, 'mmod_input.maegz')
    sdf_file = os.path.join(output_path, 'mmod_input.sdf')
    if os.path.isfile(mae_file):
        os.remove(mae_file)

    mols = mols_from_paths(input_paths, recursive=recursive)
    with open(sdf_file, 'w') as file:
        for i, mol in enumerate(mols):
            if opls_dir is None:
                mol.properties['method'] = 'opls3e_default'
                mol.properties['basis'] = 'opls'
                mol.properties['program'] = 'macromodel'
            else:
                mol.properties['method'] = 'opls3e_custom'
                mol.properties['basis'] = 'opls'
                mol.properties['program'] = 'macromodel'
            mol.to_file(file, 'SDF')
            
    sdf_to_mae(sdf_file, mae_file, schrodinger_path)
    
    com_file = os.path.join(output_path, 'mmod_openff_benchmark.com')
    with open(com_file, 'w') as outfile:
        outfile.write('mmod_input.maegz\n')
        outfile.write('mmod_output.maegz\n')
        outfile.write(' MMOD       0      1      0      0     0.0000     0.0000     0.0000     0.0000\n'\
                      ' DEBG      55      0      0      0     0.0000     0.0000     0.0000     0.0000\n'\
                      ' FFLD      16      1      0      0     1.0000     0.0000     0.0000     0.0000\n'\
                      ' BDCO       0      0      0      0    41.5692 99999.0000     0.0000     0.0000\n'\
                      ' CRMS       0      0      0      0     0.0000     0.5000     0.0000     0.0000\n'\
                      ' BGIN       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n'\
                      ' READ       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n'\
                      ' CONV       1      0      0      0     5.0e-9     0.0000     0.0000     0.0000\n'\
                      ' MINI      10      0   1500      0     0.0000     0.0000     0.0000     0.0000\n'\
                      ' END        0      0      0      0     0.0000     0.0000     0.0000     0.0000\n')

        
    command = [os.path.join(schrodinger_path, 'macromodel'),
               'mmod_openff_benchmark.com'
           ]
    if host_settings is not None:
        command.append('-HOST')
        command.append(f'{host_settings}')
    if opls_dir is not None:
        command.append('-OPLSDIR')
        command.append(os.path.abspath(opls_dir))
    subprocess.run(command, cwd=os.path.join(output_path))



def postprocess(input_paths, schrodinger_path, output_path):
    """Postprocess output from Schrodinger macromodel OPLS3e minimization.
    
    Parameters
    ----------
    input_paths : str or list of str
        Schrodinger *.mae or *.maegz file (output from schrodinger optimize).
    schrodinger_path : str
        Path to Schrodinger executables
    output_path: str
        Path where output sdf files should be written to    
    """
    
    mols = []
    for path in input_paths:
        if (os.path.isfile(path) and (path.split('.')[-1].lower() == 'maegz' or path.split('.')[-1].lower() == 'mae')):
            sdf_file = '.'.join(path.split('.')[:-1]+['sdf'])
            mae_to_sdf(path, sdf_file, schrodinger_path, append=False)
            tmp_mols = Molecule.from_file(sdf_file, 'SDF')
            if type(tmp_mols) == Molecule:
                mols.append(tmp_mols)
            else:
                mols += tmp_mols

    os.makedirs(output_path, exist_ok=True)
    os.makedirs(os.path.join(output_path, 'opls3e_custom'), exist_ok=True)
    os.makedirs(os.path.join(output_path, 'opls3e_default'), exist_ok=True)

    for mol in mols:
        mol.properties.pop('initial_energy')
        mol.properties['final_energy'] = f'{mol.properties["r_mmod_Potential_Energy-OPLS3e"]} kJ / mole'
        mol.properties['molecule_index'] = f"{mol.properties['molecule_index']:05d}"
        mol.properties['conformer_index'] = f"{mol.properties['conformer_index']:02d}"
        if mol.properties['method'] != 'opls3e_default' and mol.properties['method'] != 'opls3e_custom':
            raise ValueError(f'Method {mol.properties["method"]} not known.')
        else:
            mol.to_file(os.path.join(output_path, mol.properties['method'], f'{mol.name}.sdf'), 'SDF')

