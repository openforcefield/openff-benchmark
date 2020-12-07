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
from rdkit import Chem
from rdkit.Chem import rdMolAlign

from . import metrics, readwrite


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


def main(input_path, ref_method, output_directory="./results"):
    """
    For 2+ sets of SDF files that are analogous in terms of molecules and their
    conformers, assess them with respective to a reference SDF file (e.g., QM).
    Metrics include RMSD of conformers, TFD, and relative energy differences.

    Parameters
    ----------
    input_path : str
        Path to directory with SDF files of molecules. 
        Multiple input paths can be specified.
    ref_method : str
        Tag of reference methods. The molecules having this tag in
        the "method" SDF property will be used as reference
    output_directory : str
        Directory path to deposit exported data. If not present, this 
        directory will be created. default: ./results/
    """
    mols = {}
    dataframes = {}
    for path in input_path:
        # read in molecules
        m_mols = readwrite.read_sdfs(path)

        # specify method of molecules
        method = m_mols[0].properties['method']
        # assert that all molecules of path are from the same method
        for mol in m_mols:
            assert mol.properties["method"] == method, f"Molecules of different methods in path {path}."

        # append method, mols and dataframes
        mols[method] = m_mols
        dataframes[method] = readwrite.mols_to_dataframe(mols[method])

    assert ref_method in mols, f"Input path for reference method {ref_method} not specified."
    
    ref_confs = get_ref_confs(dataframes[ref_method])

    os.makedirs(output_directory, exist_ok=True)
    for m in dataframes:
        ref_to_ref_confs(dataframes[m], ref_confs)
        calc_rmsd(dataframes[ref_method], dataframes[m])
        calc_tfd(dataframes[ref_method], dataframes[m])
        calc_dde(dataframes[ref_method], dataframes[m])

        readwrite.write_results(dataframes[m], os.path.join(output_directory, f"{m}.csv"))

### ------------------- Parser -------------------

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    # run main
    print("Log file from compare_ffs.py")
    main()
