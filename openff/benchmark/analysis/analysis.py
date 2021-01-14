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
import warnings

from rdkit import Chem
from rdkit.Chem import rdMolAlign
from tqdm import tqdm

from . import metrics, readwrite

def get_ref_confs(dataframe):
    ref_confs = {}
    dataframe.loc[:,'ref_conf'] = False
    for mid in tqdm(dataframe.molecule_index.unique(), desc='Finding reference molecules'):
        confs = dataframe.loc[dataframe.molecule_index==mid]
        if confs.shape[0] == 1:
            ref_conf = confs.name[0]
        else:
            ref_conf = confs.final_energy.idxmin()
        ref_confs[mid] = ref_conf
        dataframe.loc[ref_conf, 'ref_conf'] = True
    return ref_confs

def ref_to_ref_confs(dataframe, ref_confs):
    for mid in tqdm(dataframe.molecule_index.unique(), desc='Referencing energies'):
        confs = dataframe.loc[dataframe.molecule_index==mid]
        ref_conf = ref_confs[mid]
        ref_energy = confs.loc[ref_conf, 'final_energy']
        for i, row in confs.iterrows():
            dataframe.loc[i, 'final_energy'] = row['final_energy'] - ref_energy
            
def calc_tfd(reference, result):
    for i, row in tqdm(reference.iterrows(), desc='Calculating TFD'):
        result.loc[i, 'tfd'] = metrics.calc_tfd(row['mol'].to_rdkit(), result.loc[i, 'mol'].to_rdkit())

def calc_rmsd(reference, result):
    for i, row in tqdm(reference.iterrows(), desc='Calculating RMSD'):
        result.loc[i, 'rmsd'] = rdMolAlign.GetBestRMS(row['mol'].to_rdkit(), result.loc[i, 'mol'].to_rdkit())
        
def calc_dde(reference, result):
    result.loc[:,'dde'] = result.final_energy - reference.final_energy

def match_minima(input_path, ref_method, output_directory="./results"):
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
    ref_to_ref_confs(dataframes[ref_method], ref_confs)
    os.makedirs(output_directory, exist_ok=True)
    matches = {}
    for m in dataframes:
        # if m == ref_method:
        #     continue
        match = compare_conformers(dataframes[ref_method], dataframes[m], 1.0)
        new_ref_confs = {i: match[ match['name'] == ref_conf ]['ff_mol_name'].values[0] for i, ref_conf in ref_confs.items() }
        ref_to_ref_confs(dataframes[m], new_ref_confs)
        for i, row in match.iterrows():
            match.loc[i, 'tfd'] = metrics.calc_tfd(dataframes[ref_method].loc[row['name'], 'mol'].to_rdkit(), dataframes[m].loc[row['ff_mol_name'], 'mol'].to_rdkit())
            match.loc[i, 'dde'] = dataframes[m].loc[row['ff_mol_name'], 'final_energy'] - dataframes[ref_method].loc[row['name'], 'final_energy']
        matches[m] = match
        readwrite.write_results(match, 
                                os.path.join(output_directory, f"matched_{m}.csv"), 
                                columns=['name', 'group_name', 'molecule_index', 'conformer_index', 'ff_mol_name', 'rmsd', 'tfd', 'dde']
        )

def compare_conformers(reference, result, rmsd_cutoff):
    """
    For different methods, match the conformer minima to those of the reference
    method. Ex. Conf G of reference method matches with conf R of method 2.

    Parameters
    ----------
    in_dict : OrderedDict
        dictionary from input file, where key is method and value is dictionary
        first entry should be reference method
        in sub-dictionary, keys are 'sdfile' and 'sdtag'
    rmsd_cutoff : float
        cutoff above which two structures are considered diff conformers

    Returns
    -------
    mol_dict : dict of dicts
        mol_dict['mol_name']['energies'] =
            [[file1_conf1_E file1_conf2_E] [file2_conf1_E file2_conf2_E]]
        An analogous structure is followed for mol_dict['mol_name']['indices'].

    """
    
    conformer_match = reference.copy()
    for mid in tqdm(reference.molecule_index.unique(), desc='Matching conformers'):
        ref_confs = reference.loc[reference.molecule_index==mid]
        query_confs = result.loc[result.molecule_index==mid]
        rms_matrix = {i: {} for i, ref_row in ref_confs.iterrows()}
        for i, ref_row in ref_confs.iterrows():
            for j, query_row in query_confs.iterrows():
                rmsd = rdMolAlign.GetBestRMS(ref_row['mol'].to_rdkit(), query_row['mol'].to_rdkit())
                rms_matrix[i][j] = rmsd
        for ref, rms_list in rms_matrix.items():
            conf = min(rms_list, key=rms_list.get)
            conformer_match.loc[ref, 'ff_mol_name'] = conf
            conformer_match.loc[ref, 'rmsd'] = rms_list[conf]

    return conformer_match



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
    for path in tqdm(input_path, desc='Reading files'):
        # read in molecules
        m_mols = readwrite.read_sdfs(path)

        # assert that all molecules of path are from the same method
        for mol in m_mols:
            # specify method of molecules
            method = mol.properties['method']
            if method in mols:
                mols[method].append(mol)
            else:
                mols[method] = [ mol ]
                
    # convert molecules to dataframes
    for method in mols:
        dataframes[method] = readwrite.mols_to_dataframe(mols[method])

    assert ref_method in mols, f"Input path for reference method {ref_method} not specified."

    index_intersect = dataframes[ref_method].index
    for m in tqdm(dataframes, desc='Checking input'):
        index_intersect = index_intersect.intersection(dataframes[m].index)
    for m, df in tqdm(dataframes.items(), desc='Checking input'):
        dataframes[m] = df.loc[index_intersect]
        if dataframes[m].shape != df.shape:
            warnings.warn(f"Not all conformers of method {m} considered, because these are not available in other methods.")


    ref_confs = get_ref_confs(dataframes[ref_method])
    ref_to_ref_confs(dataframes[ref_method], ref_confs)

    os.makedirs(output_directory, exist_ok=True)
    for m in tqdm(dataframes, desc='Processing data'):
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
