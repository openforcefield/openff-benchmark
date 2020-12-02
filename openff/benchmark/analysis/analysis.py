#!/usr/bin/env python

"""
analysis.py

For 2+ sets of SDF files that are analogous in terms of molecules and their conformers,
assess them (e.g., having FF geometries) with respective to a reference SDF
file (e.g., having QM geometries). Metrics include: RMSD of conformers, TFD
(another geometric evaluation), and relative energy differences.

By:      David F. Hahn
Version: Dec 1 2020

"""

import os
import numpy as np
import pandas as pd
import metrics
import readwrite
from rdkit import Chem
from rdkit.Chem import rdMolAlign

def get_ref_confs(dataframe):
    ref_confs = {}
    dataframe.loc[:,'ref_conf'] = False
    for mid in dataframe.molecule_index.unique():
        confs = dataframe.loc[dataframe.molecule_index==mid]
        if confs.shape[0] == 1:
            ref_conf = confs.name[0]
        else:
            ref_conf = confs.final_energy.idxmin()
        ref_confs[mid] = ref_conf
        dataframe.loc[ref_conf, 'ref_conf'] = True
    return ref_confs

def ref_to_ref_confs(dataframe, ref_confs):
    for mid in dataframe.molecule_index.unique():
        confs = dataframe.loc[dataframe.molecule_index==mid]
        ref_conf = ref_confs[mid]
        for i, row in confs.iterrows():
            dataframe.loc[i, 'final_energy'] -= dataframe.loc[ref_conf, 'final_energy']
            
def calc_tfd(reference, result):
    for i, row in reference.iterrows():
        result.loc[i, 'tfd'] = metrics.calc_tfd(row['mol'].to_rdkit(), result.loc[i, 'mol'].to_rdkit())

def calc_rmsd(reference, result):
    for i, row in reference.iterrows():
        result.loc[i, 'rmsd'] = rdMolAlign.GetBestRMS(row['mol'].to_rdkit(), result.loc[i, 'mol'].to_rdkit())
        
def calc_dde(reference, result):
    result.loc[:,'dde'] = result.final_energy - reference.final_energy


def main():
    """
    For 2+ SDF files that are analogous in terms of molecules and their
    conformers, assess them with respective to a reference SDF file (e.g., QM).
    Metrics include RMSD of conformers, TFD, and relative energy differences.

    Parameter
    ---------
    in_dict : OrderedDict
        dictionary from input file, where key is method and value is dictionary
        first entry should be reference method
        in sub-dictionary, keys are 'sdfile' and 'sdtag'
    read_pickle : Boolean
        read in data from metrics.pickle
    conf_id_tag : string
        label of the SD tag that should be the same for matching conformers
        in different files
    plot : Boolean
        generate line plots of conformer energies
    mol_slice : numpy slice object
        The resulting integers are numerically sorted and duplicates removed.
        e.g., slices = np.s_[0, 3:5, 6::3] would be parsed to return
        [0, 3, 4, 6, 9, 12, 15, 18, ...]
        Can also parse from end: [-3:] gets the last 3 molecules, and
        [-2:-1] is the same as [-2] to get just next to last molecule.

    """
    # for now hardcoded, could be taken from SEASONS of David D. or from CL arg
    methods = ['qcarchive', 'openff-1.2.0']
    ref_method = 'qcarchive'

    mols = {}
    dataframes = {}
    for m in methods:
        mols[m] = readwrite.read_sdfs('qcarchive')
        dataframes[m] = readwrite.mols_to_dataframe(mols[m])

    ref_confs = get_ref_confs(dataframes[ref_method])

    for m in methods:
        ref_to_ref_confs(dataframes[m], ref_confs)
        if m == ref_method:
            continue
        else:
            calc_rmsd(dataframes[ref_method], dataframes[m])
            calc_tfd(dataframes[ref_method], dataframes[m])
            calc_dde(dataframes[ref_method], dataframes[m])

    os.makedirs('results', exist_ok=True)
    for m in methods:
        if m == ref_method:
            continue
        else:
            dataframes[m].to_csv(f'results/{m}.csv')

### ------------------- Parser -------------------

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    # run main
    print("Log file from compare_ffs.py")
    main()
