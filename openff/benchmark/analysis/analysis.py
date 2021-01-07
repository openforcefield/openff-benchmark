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

    os.makedirs(output_directory, exist_ok=True)
    for m in dataframes:
        if m == ref_method:
            continue
        match = compare_conformers(dataframes[ref_method], dataframes[m], 1.0)
        new_ref_confs = {i: match[ match['Ref_mol']== ref_conf ]['FF_mol'].values[0] for i, ref_conf in ref_confs.items() }
        ref_to_ref_confs(dataframes[m], new_ref_confs)
        for i, row in match.iterrows():
            match.loc[i, 'tfd'] = metrics.calc_tfd(dataframes[ref_method].loc[row['Ref_mol'], 'mol'].to_rdkit(), dataframes[m].loc[row['FF_mol'], 'mol'].to_rdkit())
            match.loc[i, 'dde'] = dataframes[m].loc[row['FF_mol'], 'final_energy'] - dataframes[ref_method].loc[row['Ref_mol'], 'final_energy']

        match.to_csv(os.path.join(output_directory, f"match_{m}.csv"),
                     index=False, 
                     float_format='%15.8e'
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
    
    conformer_match=[]
    for mid in reference.molecule_index.unique():
        ref_confs = reference.loc[reference.molecule_index==mid]
        query_confs = result.loc[result.molecule_index==mid]
        rms_matrix = {i: {} for i, ref_row in ref_confs.iterrows()}
        for i, ref_row in ref_confs.iterrows():
            for j, query_row in query_confs.iterrows():
                rmsd = rdMolAlign.GetBestRMS(ref_row['mol'].to_rdkit(), query_row['mol'].to_rdkit())
                rms_matrix[i][j] = rmsd
        for ref, rms_list in rms_matrix.items():
            conf = min(rms_list, key=rms_list.get)
            conformer_match.append([ref, conf, rms_list[conf]])

    match = pd.DataFrame(conformer_match, columns=['Ref_mol', 'FF_mol', 'rmsd'])

    return match


    # assess each file against reference
    for ff_label, ff_dict in in_dict.items():
        sdf_query = ff_dict['sdfile']
        sdf_tag = ff_dict['sdtag']

        # load molecules from open reference and query files
        print(f"\n\nOpening reference file {sdf_ref}")
        mols_ref = reader.read_mols(sdf_ref)

        print(f"Opening query file {sdf_query} for [ {ff_label} ] energies")
        mols_query = reader.read_mols(sdf_query)

        # loop over each molecule in reference and query files
        for rmol in mols_ref:
            mol_name = rmol.GetTitle()
            ref_nconfs = rmol.NumConfs()
            run_match = False

            for qmol in mols_query:

                # same mol titles should mean same molecular identity;
                # when same molecular identity found, break out of loop to
                # start matching conformers
                if rmol.GetTitle() == qmol.GetTitle():
                    run_match = True
                    break

            # create entry for this mol in mol_dict if not already present
            # energies [i][j] will be 2d list of ith method and jth conformer
            if mol_name not in mol_dict:
                mol_dict[mol_name] = {'energies': [], 'indices': []}

            # no same molecules were found bt ref and query methods
            # for N reference minima of each mol, P matching indices for each ref minimia
            if not run_match:
                print(f"No \"{mol_name}\" molecule found in {sdf_query}")

                # fill in -2 error values for conformer indices
                mol_dict[mol_name]['indices'].append([-2] * ref_nconfs)

                # fill in nan values for conformer energies and ref_nconfs
                mol_dict[mol_name]['energies'].append([np.nan] * ref_nconfs)

                # reset mols_query generator
                mols_query = reader.read_mols(sdf_query)

                # continue with the next rmol
                continue

            # get data from specified sd tag for all conformers
            data_confs = reader.get_sd_list(qmol, sdf_tag)

            # format sd tag data to float types
            float_data_confs = list(map(float, data_confs))

            # store sd data from tags into dictionary
            mol_dict[mol_name]['energies'].append(float_data_confs)

            # don't run match if query method is same as reference method
            # keep this section after sd tag extraction of energies
            if sdf_query == sdf_ref:
                print("Skipping comparison against self.")
                mol_dict[mol_name]['indices'].append([-1] * ref_nconfs)
                continue

            # run the match here
            # get indices of qmol conformers that match rmol conformers
            molIndices = compare_two_mols(rmol, qmol, rmsd_cutoff)
            mol_dict[mol_name]['indices'].append(molIndices)

    return mol_dict

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
