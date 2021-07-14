#!/usr/bin/env python

"""
new_analysis.py

For 2+ sets of SDF files that are analogous in terms of molecules and their conformers,
assess them (e.g., having FF geometries) with respective to a reference SDF
file (e.g., having QM geometries). Metrics include: RMSD of conformers, and relative energy differences.

By:      David F. Hahn, Lorenzo D'Amore
Version: May 6 2021

"""

import os
import numpy as np
import pandas as pd
import warnings

from rdkit import Chem
from rdkit.Chem import rdMolAlign
from tqdm import tqdm

from . import metrics, readwrite

def lucas(input_path, ref_method, output_directory="./5-results-lucas"):
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

    dataframe = dataframes[ref_method]

    # find QM minimas qm_mins
    qm_mins = {}
    dataframe.loc[:,'qm_min_check'] = False
    for mid in tqdm(dataframe.molecule_index.unique(), desc='Finding QM minima'):
        confs = dataframe.loc[dataframe.molecule_index==mid]
        if confs.shape[0] == 1:
            qm_min = confs.name[0]
        else:
            qm_min = confs.final_energy.idxmin()
        qm_mins[mid] = qm_min
        dataframe.loc[qm_min, 'qm_min_check'] = True

    os.makedirs(output_directory, exist_ok=True)
    matches = {}

    # find the MM conformer closest to QM qm_min based on rmsd
    # loop over methods

    for m in dataframes:
        
        # if the method is the reference method, we do not do the comparison
        # because it's just a comparison with itself
        if m == ref_method:
            continue

        mm_dataframe = dataframes[m]    
                                                                                                            
        match = compare_conformers(dataframes[ref_method], dataframes[m])
        mm_ref_confs = {molecule_id: 
                                    match[ match['name'] == ref_conformer ]['ff_mol_name'].values[0]
                                    for molecule_id, ref_conformer in qm_mins.items() }
        matches[m] = match
                                                                                                                                                                            
        # find MM minimas mm_mins
        mm_mins = {}
        mm_dataframe.loc[:,'mm_min_check'] = False
        for mid in tqdm(mm_dataframe.molecule_index.unique(), desc='Finding QM minima'):
            confs = mm_dataframe.loc[mm_dataframe.molecule_index==mid]
            if confs.shape[0] == 1:
                mm_min = confs.name[0]
            else:
                mm_min = confs.final_energy.idxmin()
            mm_mins[mid] = mm_min
            mm_dataframe.loc[mm_min, 'mm_min_check'] = True

        # calculate dE between mm_ref_conf and mm_min
        for i, row in mm_dataframe.iterrows():
            mm_min = mm_mins[row['molecule_index']]
            mm_ref_conf = mm_ref_confs[row['molecule_index']]
            if row['mm_min_check'] == True:
                mm_dataframe.loc[i, 'final_energy'] = mm_dataframe.loc[mm_ref_conf,'final_energy'] - mm_dataframe.loc[mm_min,'final_energy']

        # calculate RMSD between mm_min and mm_ref_conf
        for i, row in mm_dataframe.iterrows():
            mm_min = mm_mins[row['molecule_index']]
            mm_ref_conf = mm_ref_confs[row['molecule_index']]
            if row['mm_min_check'] == True:
                mm_dataframe.loc[i, 'rmsd'] = rdMolAlign.GetBestRMS(
                                dataframe.loc[mm_ref_conf, 'mol'].to_rdkit(), dataframe.loc[mm_min, 'mol'].to_rdkit())
        
        # take only the mm_min = True of the dataframe
        mm_results = mm_dataframe.loc[mm_dataframe.mm_min_check]
        
        # adds qm_min and mm_ref_conf to the new dataframe
        for i, row in mm_results.iterrows():
            qm_min = qm_mins[row['molecule_index']]
            mm_ref_conf = mm_ref_confs[row['molecule_index']]
            mm_results.loc[i, 'qm_min'] = qm_min
            mm_results.loc[i, 'mm_ref_conf'] = mm_ref_conf

        mm_results = mm_results.rename(columns={'name':'mm_min', 'mm_ref_conf':'mm_ref', 
                                                'rmsd':'rmsd (mm_ref / mm_min)', 'final_energy':'dE (mm_ref - mm_min)'})
        
        mm_results = mm_results[['qm_min', 'mm_ref', 'mm_min', 'rmsd (mm_ref / mm_min)', 'dE (mm_ref - mm_min)']]
        
        mm_results.to_csv(os.path.join(output_directory, f"lucas_{m}.csv"), index=False, float_format='%15.8e')


def compare_conformers(reference, result):
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


def swope(input_path, ref_method, output_directory="./5-results-swope"):
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

    dataframe = dataframes[ref_method]
    
    # find QM minimas qm_mins
    qm_mins = {}
    dataframe.loc[:,'qm_min_check'] = False
    for mid in tqdm(dataframe.molecule_index.unique(), desc='Finding QM minima'):
        confs = dataframe.loc[dataframe.molecule_index==mid]
        if confs.shape[0] == 1:
            qm_min = confs.name[0]
        else:
            qm_min = confs.final_energy.idxmin()
        qm_mins[mid] = qm_min
        dataframe.loc[qm_min, 'qm_min_check'] = True

    os.makedirs(output_directory, exist_ok=True)

    # loop over methods
    for m in dataframes:
        # if the method is the reference method, we do not do the comparison
        # because it's just a comparison with itself
        if m == ref_method:
            continue

        mm_dataframe = dataframes[m]

        # find MM minimas mm_mins
        mm_mins = {}
        mm_dataframe.loc[:,'mm_min_check'] = False
        for mid in tqdm(mm_dataframe.molecule_index.unique(), desc='Finding MM minima'):
            confs = mm_dataframe.loc[mm_dataframe.molecule_index==mid]
            if confs.shape[0] == 1:
                mm_min = confs.name[0]
            else:
                mm_min = confs.final_energy.idxmin()
            mm_mins[mid] = mm_min
            mm_dataframe.loc[mm_min, 'mm_min_check'] = True

        # MM confs dE with respect to mm_min conf        
        for mid in tqdm(mm_dataframe.molecule_index.unique(), desc='Referencing MM energies'):
            confs = mm_dataframe.loc[mm_dataframe.molecule_index==mid]
            mm_min = mm_mins[mid]
            mm_min_energy = confs.loc[mm_min, 'final_energy']
            for i, row in confs.iterrows():
                mm_dataframe.loc[i, 'final_energy'] = row['final_energy'] - mm_min_energy

        # calculate RMSD
        for i, row in tqdm(mm_dataframe.iterrows(), desc='Calculating RMSD'):
            qm_min = qm_mins[row['molecule_index']]
            mm_dataframe.loc[i, 'rmsd'] = rdMolAlign.GetBestRMS(
                             row['mol'].to_rdkit(), dataframe.loc[qm_min, 'mol'].to_rdkit())

        mm_results = mm_dataframe.copy()

        # adds qm_min and mm_minf to the new dataframe
        for i, row in mm_results.iterrows():
            qm_min = qm_mins[row['molecule_index']]
            mm_min = mm_mins[row['molecule_index']]
            mm_results.loc[i, 'qm_min'] = qm_min
            mm_results.loc[i, 'mm_min'] = mm_min

        mm_results = mm_results.rename(columns={'name':'mm_conf', 'rmsd':'rmsd (qm_min / mm_conf)', 
                                                'final_energy':'dE (mm_conf - mm_min)'})
        mm_results = mm_results[['qm_min', 'mm_conf', 'mm_min', 'rmsd (qm_min / mm_conf)', 
                                'dE (mm_conf - mm_min)']]

        mm_results.to_csv(os.path.join(output_directory, f"swope_{m}.csv"), index=False, float_format='%15.8e')

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
