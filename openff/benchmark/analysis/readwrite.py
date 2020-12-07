#!/usr/bin/env python

"""
io.py

I/O operations for the analysis/report part of the openff-benchmark workflow

By:      David F. Hahn
Version: Nov 18 2020

"""

import os
import numpy as np
import pandas as pd
    
import pint
ureg = pint.UnitRegistry()

from openforcefield.topology import Molecule
from openforcefield.utils.toolkits import GLOBAL_TOOLKIT_REGISTRY, OpenEyeToolkitWrapper


oetk_loaded = False
for tkw in GLOBAL_TOOLKIT_REGISTRY.registered_toolkits:
    if isinstance(tkw, OpenEyeToolkitWrapper):
        oetk_loaded = True
if oetk_loaded:
    GLOBAL_TOOLKIT_REGISTRY.deregister_toolkit(OpenEyeToolkitWrapper)

def read_sdfs(path):
    mols = []
    if os.path.isdir(path):
        for root, dirs, files in os.walk(path):
            for file in files:
                file_name = os.path.join(root, file)
                if (os.path.exists(file_name) and file_name.split('.')[-1] == 'sdf' ):
                    mols.append(Molecule.from_file(file_name, 'SDF', allow_undefined_stereo=True))
    return mols

def convert_to_quantity(dataframe, columns='final_energy', to='kilocalories / mole'):
    if type(columns) is str:
        columns=[columns]
    for col in columns:
        dataframe[col] = dataframe[col].apply(lambda val: ureg.Quantity(val).to(to).magnitude)

def mols_to_dataframe(mols):
    moldata = []
    for mol in mols:
        moldata_i = {
                'name': mol.name,
                'mol': mol
                }
        for key, item in mol.properties.items():
            moldata_i[key] = item

        moldata.append(moldata_i)

    df = pd.DataFrame(moldata)
    df['molecule_index'] = df['molecule_index'].astype(int)
    df['conformer_index'] = df['conformer_index'].astype(int)
    convert_to_quantity(df, columns=['initial_energy', 'final_energy'])
    df.set_index('name', drop=False, inplace=True)
    return df


def write_results(dataframe, file_name, columns=['name', 'group_name', 'molecule_index', 'conformer_index', 'rmsd', 'tfd', 'dde']):
    if 'molecule_index' in columns:
        dataframe['molecule_index'] = dataframe['molecule_index'].apply(lambda x: f'{x:05d}')
    if 'conformer_index' in columns:
        dataframe['conformer_index'] = dataframe['conformer_index'].apply(lambda x: f'{x:02d}')
    dataframe[columns].to_csv(file_name, 
                              index=False if 'name' in columns else True, 
                              float_format='%15.8e')
    
