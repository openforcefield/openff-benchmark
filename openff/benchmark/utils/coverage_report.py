"""
Functions to help with producing simple coverage reports for a set of molecules.
"""
from typing import Dict, Union, List, Optional
from openforcefield.typing.engines.smirnoff import ForceField
from openforcefield.topology import Molecule
import glob
import os
import tqdm
import shutil
import logging

logger = logging.getLogger('openforcefield.utils.toolkits')
prev_log_level = logger.getEffectiveLevel()
logger.setLevel(logging.ERROR)

from openforcefield.topology import Molecule
from openforcefield.utils.toolkits import GLOBAL_TOOLKIT_REGISTRY, OpenEyeToolkitWrapper

logger.setLevel(prev_log_level)

oetk_loaded = False
for tkw in GLOBAL_TOOLKIT_REGISTRY.registered_toolkits:
    if isinstance(tkw, OpenEyeToolkitWrapper):
        oetk_loaded = True
if oetk_loaded:
    GLOBAL_TOOLKIT_REGISTRY.deregister_toolkit(OpenEyeToolkitWrapper)

def _update_coverage(main_report, single_report):
    """
    Take the single molecule report and update the master coverage dict
    """
    for param_type in main_report.keys():
        for param_id, no in single_report[param_type].items():
            if param_id not in main_report[param_type]:
                # set the param number if new
                main_report[param_type][param_id] = no
            else:
                # add onto the current parameter count
                main_report[param_type][param_id] += no

def generate_coverage_report(input_molecules: List[Molecule],
                             forcefield_name: str,
                             processors: Optional[int] = None):
    """
    For the given set of molecules produce a coverage report of the smirks parameters used in typing the molecule.
    Also try and produce charges for the molecules as some may have missing bccs.

    Parameters
    ----------
    input_molecules: The List of 3d sdf files to run the coverage on or the name of the directory containing the files.
    forcefield_name: The name of the openFF forcefield to run the coverage for.
    processors: The number of processors we can use to build the coverage report

    Returns
    -------
    coverage_report: A dictionary split into parameter types which lists the number of occurrences of each parameter.
    success_mols: A list of openforcefield.topology.Molecule objects that were successful in this step
    error_mols: A list of tuples (Molecule, Exception) of molecules that failed this step
    """
    if isinstance(input_molecules, Molecule):
        input_molecules = [input_molecules, ]
    coverage = {"Angles": {}, "Bonds": {}, "ProperTorsions": {}, "ImproperTorsions": {}, "vdW": {}}


    # make the forcefield
    ff = ForceField(forcefield_name)
    # For speed, don't test charge assignment for now
    ff.deregister_parameter_handler('ToolkitAM1BCC')
    ff.get_parameter_handler('ChargeIncrementModel',
                             {'partial_charge_method':'formal_charge',
                              'version':'0.3'})
    # now run coverage on each molecule
    success_mols = []
    error_mols = []
    if processors is None or processors > 1:
        from multiprocessing import Pool

        with Pool(processes=processors) as pool:
            # generate the work list
            work_list = [pool.apply_async(single_molecule_coverage, (molecule, ff)) for molecule in input_molecules]
            for work in tqdm.tqdm(
                    work_list,
                    total=len(work_list),
                    ncols=80,
                    desc="{:30s}".format("Generating Coverage"),
            ):
                report, e = work.get()
                molecule = report["molecule"]
                if e is None:
                    _update_coverage(coverage, report)
                    success_mols.append(molecule)
                #except Exception as e:
                else:
                    error_mols.append((molecule, e))

    else:
        for molecule in tqdm.tqdm(input_molecules,
                                  total=len(input_molecules),
                                  ncols=80,
                                  desc="{:30s}".format("Generating Coverage")):
            #try:
            report, e = single_molecule_coverage(molecule, ff)
            # update the master coverage dict
            mol: Molecule = report["molecule"]
            if e is None:
                _update_coverage(coverage, report)
                success_mols.append(mol)
            #except Exception as e:
            else:
                error_mols.append((molecule, e))

    # Sort the keys of the coverage dict, so that it doesn't show which parameters were found first
    for parameter_type in coverage.keys():
        sorted_sub_dict = dict(sorted(coverage[parameter_type].items()))
        coverage[parameter_type] = sorted_sub_dict

    # record how many molecules we processed
    coverage["passed_unique_molecules"] = len(success_mols)
    coverage["total_unique_molecules"] = len(success_mols) + len(error_mols)
    coverage["forcefield_name"] = forcefield_name

    return coverage, success_mols, error_mols


def single_molecule_coverage(molecule: Molecule,
                             forcefield: ForceField):
    #-> Dict[str, Dict[str, int]], List[Molecule], List[Molecule, Exception]
    """
    For a single molecule generate a coverage report and try to build an openmm system this will also highlight any missing parameters and dificulties with charging the molecule.

    Parameters
    ----------
    molecule: The openforcefield molecule object for which the report should be generated.
    ff: The openforcefield typing engine that should be used to check coverage and build an openmm system.

    Returns
    -------
    report: dict
        A dictionary of the coverage report 
    e: Exception or None. 
        The exception raised in this step, if any. 
        If not None, it should be assumed that coverage is invalid.
    """

    coverage = {"Angles": {}, "Bonds": {}, "ProperTorsions": {}, "ImproperTorsions": {}, "vdW": {}}
    coverage["molecule"] = molecule
    try:
        labels = forcefield.label_molecules(molecule.to_topology())[0]
        for param_type, params in labels.items():
            for param in params.values():
                if param.id not in coverage[param_type]:
                    coverage[param_type][param.id] = 1
                else:
                    coverage[param_type][param.id] += 1
        # now generate a system this will catch any missing parameters
        # and molecules that can not be charged
        _ = forcefield.create_openmm_system(molecule.to_topology())
        return coverage, None
    except Exception as e:
        return coverage, e

