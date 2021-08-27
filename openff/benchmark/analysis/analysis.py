#!/usr/bin/env python

"""
analysis.py

For 2+ sets of SDF files that are analogous in terms of molecules and their conformers,
assess them (e.g., having FF geometries) with respective to a reference SDF
file (e.g., having QM geometries). Metrics include: RMSD of conformers, TFD
(another geometric evaluation), and relative energy differences.

By:      David F. Hahn, Lorenzo D'Amore
Version: Jul 22 2021

"""

import os
import numpy as np
import pandas as pd
import warnings

from rdkit import Chem
from rdkit.Chem import rdMolAlign
from tqdm import tqdm

from . import metrics, readwrite


def gather_df(input_path, ref_method):
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

    assert ref_method in mols, f"No molecules for reference method {ref_method} in input path(s)."

    return dataframes


def intersect(dataframes, ref_method):
    index_intersect = dataframes[ref_method].index
    for m in tqdm(dataframes, desc='Checking input'):
        index_intersect = index_intersect.intersection(dataframes[m].index)
    for m, df in tqdm(dataframes.items(), desc='Checking input'):
        dataframes[m] = df.loc[index_intersect]
        if dataframes[m].shape != df.shape:
            warnings.warn(f"Not all conformers of method {m} considered, because these are not available in other methods.")


def get_confs_min(dataframe):
    confs_min = {}
    dataframe.loc[:,'conf_min'] = False
    for mid in tqdm(dataframe.molecule_index.unique(), desc='Finding conformer minima'):
        confs = dataframe.loc[dataframe.molecule_index==mid]
        if confs.shape[0] == 1:
            conf_min = confs.name[0]
        else:
            conf_min = confs.final_energy.idxmin()
        confs_min[mid] = conf_min
        dataframe.loc[conf_min, 'conf_min'] = True
    return confs_min


def calc_de(dataframe, confs_min):
    for mid in tqdm(dataframe.molecule_index.unique(), desc='Calculating energy difference'):
        confs = dataframe.loc[dataframe.molecule_index==mid]
        conf_min = confs_min[mid]
        ref_energy = confs.loc[conf_min, 'final_energy']
        for i, row in confs.iterrows():
            dataframe.loc[i, 'final_energy'] = row['final_energy'] - ref_energy


def calc_rmsd(reference, result, ref_name, result_name):
    for i, row in tqdm(reference.iterrows(), desc='Calculating RMSD'):
        try:
            result.loc[i, 'rmsd'] = rdMolAlign.GetBestRMS(row['mol'].to_rdkit(), result.loc[i, 'mol'].to_rdkit())
        except RuntimeError:
            result.loc[i, 'rmsd'] = np.NaN
            print(f"Unable to calculate best RMSD between {ref_name} and {result_name}; conformer `{i}`")
           

def calc_tfd(reference, result):
    for i, row in tqdm(reference.iterrows(), desc='Calculating TFD'):
        result.loc[i, 'tfd'] = metrics.calc_tfd(row['mol'].to_rdkit(), result.loc[i, 'mol'].to_rdkit())


def calc_dde(reference, result):
    result.loc[:,'dde'] = result.final_energy - reference.final_energy


def match_minima(input_path, ref_method, output_directory="./results"):
    
    # collects the input molecules and build a dataframe 
    dataframes = gather_df(input_path, ref_method)
 
    # takes only conformers which are present in all the methods
    intersect(dataframes, ref_method)

    # get a dictionary with molecule index as key and the name of the reference as item
    confs_min = get_confs_min(dataframes[ref_method])

    # reference all final energies to the reference conformer's final energy (the reference conformers final energy will be 0 afterwards)
    calc_de(dataframes[ref_method], confs_min)
    os.makedirs(output_directory, exist_ok=True)
    
    # loop over methods
    for m in dataframes:
        # if the method is the reference method, we do not do the comparison
        # because it's just a comparison with itself
        if m == ref_method:
            continue
        match = get_ref_conf(dataframes[ref_method], dataframes[m], ref_method, m)
        ref_confs = {molecule_id: 
                         match[ match['name'] == ref_conformer ]['ff_mol_name'].values[0] 
                         for molecule_id, ref_conformer in confs_min.items() }
        calc_de(dataframes[m], ref_confs)
        for i, row in match.iterrows():
            match.loc[i, 'tfd'] = metrics.calc_tfd(
                dataframes[ref_method].loc[row['name'], 'mol'].to_rdkit(), 
                dataframes[m].loc[row['ff_mol_name'], 'mol'].to_rdkit()
            )
            match.loc[i, 'dde'] = dataframes[m].loc[row['ff_mol_name'], 'final_energy'] - dataframes[ref_method].loc[row['name'], 'final_energy']
        readwrite.write_results(match, 
                                os.path.join(output_directory, f"matched_{m}.csv"), 
                                columns=['name', 'group_name', 'molecule_index', 'conformer_index', 'ff_mol_name', 'rmsd', 'tfd', 'dde']
        )


def lucas(input_path, ref_method, output_directory="./5-results-lucas"):
    """Execute comparison analysis proposed by Xavier Lucas.

    The command accepts the paths of the optimized molecules obtained from the optimization step
    and creates one output csv file per method.

    For each molecule, the code finds the MM reference conformer (ref_conf) with the lowest RMSD
    value with respect to the QM global minimum (qm_min) and then reports the relative energy (dE) 
    and RMDS between ref_conf and the MM global minimum (mm_min).

    Parameters
    ----------
    input_path : Iterable[Path-like]
        Input paths to gather input SDFs of molecule conformers to compare.
    ref_methd : str
        The value of the SDF property `method` to use as the reference method.
    output_directory : Path-like
        The directory in which to output results.

    """

    # collects the input molecules and build a dataframe
    dataframes = gather_df(input_path, ref_method)

    # takes only conformers which are present in all the methods
    intersect(dataframes, ref_method)

    # find QM min
    qm_df = dataframes[ref_method]
    qm_mins = get_confs_min(qm_df)

    os.makedirs(output_directory, exist_ok=True)

    # find the MM conformer closest to QM qm_min based on rmsd
    # loop over methods

    for m in dataframes:

        # if the method is the reference method, we do not do the comparison
        # because it's just a comparison with itself
        if m == ref_method:
            continue

        mm_df = dataframes[m]

        match = get_ref_conf(qm_df, mm_df, ref_method, m)
        ref_confs = {molecule_id:
                                 match[ match['name'] == ref_conformer ]['ff_mol_name'].values[0]
                                 for molecule_id, ref_conformer in qm_mins.items() }

        # find MM min
        mm_mins = get_confs_min(dataframes[m])

        # calculate dE between ref_conf and mm_min
        for i, row in mm_df.iterrows():
            mm_min = mm_mins[row['molecule_index']]
            ref_conf = ref_confs[row['molecule_index']]
            if row['conf_min'] == True:
                mm_df.loc[i, 'final_energy'] = mm_df.loc[ref_conf,'final_energy'] - mm_df.loc[mm_min,'final_energy']

        # calculate RMSD between mm_min and ref_conf
        for i, row in mm_df.iterrows():
            mm_min = mm_mins[row['molecule_index']]
            ref_conf = ref_confs[row['molecule_index']]
            if row['conf_min'] == True:
                try:
                    mm_df.loc[i, 'rmsd'] = rdMolAlign.GetBestRMS(
                                    qm_df.loc[ref_conf, 'mol'].to_rdkit(), qm_df.loc[mm_min, 'mol'].to_rdkit())
                except RuntimeError:
                    mm_df.loc[i, 'rmsd'] = np.NaN
                    print(f"Unable to calculate best RMSD between {ref_method} and {m}; conformer `{i}`")

        # take only the mm_min = True of the dataframe
        mm_results = mm_df.loc[mm_df.conf_min].copy()

        # adds qm_min and ref_conf to the new dataframe
        for i, row in mm_results.iterrows():
            qm_min = qm_mins[row['molecule_index']]
            ref_conf = ref_confs[row['molecule_index']]
            mm_results.loc[i, 'qm_min'] = qm_min
            mm_results.loc[i, 'ref_conf'] = ref_conf

        mm_results = mm_results.rename(columns={'name':'mm_min', 'ref_conf':'mm_ref',
                                                'rmsd':'rmsd (mm_ref/mm_min)', 'final_energy':'dE (mm_ref-mm_min)'})

        mm_results = mm_results[['qm_min', 'mm_ref', 'mm_min', 'rmsd (mm_ref/mm_min)', 'dE (mm_ref-mm_min)']]

        mm_results.to_csv(os.path.join(output_directory, f"lucas_{m}.csv"), index=False, float_format='%15.8e')


def swope(input_path, ref_method, output_directory="./5-results-swope"):
    """Execute comparison analysis proposed by William Swope.

    The command accepts the paths of the optimized molecules obtained from the optimization step
    and creates one output csv file per method.

    For each molecule, the code reports (i) the relative energy (dE) between each MM conformer and the MM 
    conformer which is the global minimum (mm_min); (ii) the RMSD between each MM conformer and the QM
    conformer which is the global minimum (qm_min).

    Parameters
    ----------
    input_path : Iterable[Path-like]
        Input paths to gather input SDFs of molecule conformers to compare.
    ref_methd : str
        The value of the SDF property `method` to use as the reference method.
    output_directory : Path-like
        The directory in which to output results.

    """

    # collects the input molecules and build a dataframe
    dataframes = gather_df(input_path, ref_method)

    # takes only conformers which are present in all the methods
    intersect(dataframes, ref_method)

    # find QM min
    qm_df = dataframes[ref_method]
    qm_mins = get_confs_min(qm_df)

    os.makedirs(output_directory, exist_ok=True)

    # loop over methods
    for m in dataframes:
        # if the method is the reference method, we do not do the comparison
        # because it's just a comparison with itself
        if m == ref_method:
            continue

        # find MM min
        mm_df = dataframes[m]
        mm_mins = get_confs_min(dataframes[m])

        # calculate dE between ech MM conformer and mm_min
        for mid in tqdm(mm_df.molecule_index.unique(), desc='Calculating energy difference'):
            confs = mm_df.loc[mm_df.molecule_index==mid]
            mm_min = mm_mins[mid]
            mm_min_energy = confs.loc[mm_min, 'final_energy']
            for i, row in confs.iterrows():
                mm_df.loc[i, 'final_energy'] = row['final_energy'] - mm_min_energy

        # calculate RMSD
        for i, row in tqdm(mm_df.iterrows(), desc='Calculating RMSD'):
            qm_min = qm_mins[row['molecule_index']]
            try:

                mm_df.loc[i, 'rmsd'] = rdMolAlign.GetBestRMS(
                                 row['mol'].to_rdkit(), qm_df.loc[qm_min, 'mol'].to_rdkit())
            except RuntimeError:
                mm_df.loc[i, 'rmsd'] = np.NaN
                print(f"Unable to calculate best RMSD between {ref_method} and {m}; conformer `{i}`")

        mm_results = mm_df.copy()

        # adds qm_min and mm_min to the new dataframe
        for i, row in mm_results.iterrows():
            qm_min = qm_mins[row['molecule_index']]
            mm_min = mm_mins[row['molecule_index']]
            mm_results.loc[i, 'qm_min'] = qm_min
            mm_results.loc[i, 'mm_min'] = mm_min

        mm_results = mm_results.rename(columns={'name':'mm_conf', 'rmsd':'rmsd (mm_conf/qm_min)',
                                                'final_energy':'dE (mm_conf-mm_min)'})
        
        mm_results = mm_results[['qm_min', 'mm_conf', 'mm_min', 'rmsd (mm_conf/qm_min)', 'dE (mm_conf-mm_min)']]

        mm_results.to_csv(os.path.join(output_directory, f"swope_{m}.csv"), index=False, float_format='%15.8e')



def get_ref_conf(reference, result, ref_name, result_name):
    """
    For each MM method, get the conformers that are the closest (by RMSD) to the global 
    minima conformers calculated with the reference (QM) method.

    Parameters
    ----------
    in_dict : OrderedDict
        dictionary from input file, where key is method and value is dictionary
        first entry should be reference method
        in sub-dictionary, keys are 'sdfile' and 'sdtag'

    Returns
    -------
    mol_dict : dict of dicts
        mol_dict['mol_name']['energies'] =
            [[file1_conf1_E file1_conf2_E] [file2_conf1_E file2_conf2_E]]
        An analogous structure is followed for mol_dict['mol_name']['indices'].

    """
    
    conformer_match = reference.copy()
    for mid in tqdm(reference.molecule_index.unique(), desc='Matching conformers'):
        confs_min = reference.loc[reference.molecule_index==mid]
        query_confs = result.loc[result.molecule_index==mid]
        rms_matrix = {i: {} for i, ref_row in confs_min.iterrows()}
        for i, ref_row in confs_min.iterrows():
            for j, query_row in query_confs.iterrows():
                try:
                    rmsd = rdMolAlign.GetBestRMS(ref_row['mol'].to_rdkit(), query_row['mol'].to_rdkit())
                except:
                    rmsd = np.NaN
                    print(f"Unable to calculate best RMSD between {ref_name} and {result_name}; conformer `{i}`")
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
    # collects the input molecules and build a dataframe
    dataframes = gather_df(input_path, ref_method)

    # takes only conformers which are present in all the methods
    intersect(dataframes, ref_method)

    confs_min = get_confs_min(dataframes[ref_method])
    calc_de(dataframes[ref_method], confs_min)

    os.makedirs(output_directory, exist_ok=True)
    for m in tqdm(dataframes, desc='Processing data'):
        if m == ref_method:
            continue
        calc_de(dataframes[m], confs_min)
        calc_rmsd(dataframes[ref_method], dataframes[m], ref_name=ref_method, result_name=m)
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
