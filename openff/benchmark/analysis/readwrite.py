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
from openforcefield.topology import Molecule


def read_sdfs(path):
    mols = []
    if os.path.isdir(path):
        for root, dirs, files in os.walk(path):
            for file in files:
                file_name = os.path.join(root, file)
                if (os.path.exists(file_name) and file_name.split('.')[-1] == 'sdf' ):
                    mols.append(Molecule.from_file(file_name, 'SDF', allow_undefined_stereo=True))
    return mols

def mols_to_dataframe(mols):
    df = pd.DataFrame()
    for i, mol in enumerate(mols):
        df.loc[i, 'name'] = mol.name
        df.loc[i, 'mol'] = mol
        for key, item in mol.properties.items():
            df.loc[i, key] = item
    df['molecule_index'] = df['molecule_index'].astype(int)
    df['conformer_index'] = df['conformer_index'].astype(int)
    df.set_index('name')
    return df
