#!/usr/bin/env python

"""
metrics.py

Metrics calculation for the analysis/report part of the openff-benchmark workflow

By:      David F. Hahn
Version: Nov 25 2020

"""

import numpy as np
from rdkit import Chem
from rdkit.Chem import TorsionFingerprints

def calc_tfd(ref_mol, query_mol):
    """
    Calculate Torsion Fingerprint Deviation between two molecular structures.
    RDKit is required for TFD calculation.

    References
    ----------
    Modified from the following code:
    https://github.com/MobleyLab/benchmarkff/03_analysis/compare_ffs.py

    TFD reference:
    https://pubs.acs.org/doi/10.1021/ci2002318

    Parameters
    ----------
    ref_mol : RDKit RDMol
    query_mol : RDKit RDMol

    Returns
    -------
    tfd : float
        Torsion Fingerprint Deviation between ref and query molecules

    """
    # check if the molecules are the same
    # tfd requires the two molecules must be instances of the same molecule
    rsmiles = Chem.MolToSmiles(ref_mol)
    qsmiles = Chem.MolToSmiles(query_mol)
    if rsmiles != qsmiles:
        print(
            f"- WARNING: The reference mol {ref_mol.GetProp('_Name')} and "
            f"query mol {query_mol.GetProp('_Name')} do NOT have the same "
            f"SMILES strings as determined by RDKit MolToSmiles. "
            f"\n {rsmiles}\n {qsmiles}"
        )
        tfd = np.nan

    # calculate the TFD
    else:
        try:
            tfd = TorsionFingerprints.GetTFDBetweenMolecules(ref_mol, query_mol)
        # triggered for molecules such as urea
        except IndexError:
            print(
                f"- Error calculating TFD on molecule {ref_mol.GetProp('_Name')}."
                " Possibly no non-terminal rotatable bonds found."
            )
            tfd = np.nan

    return tfd
