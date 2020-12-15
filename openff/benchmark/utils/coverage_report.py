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


def generate_coverage_report(input_molecules: Union[List[str], str], forcefield_name, output_directory: str = "2-coverage-report", processors: Optional[int] = None, delete_existing: bool = False) -> Dict[str, Dict[str, int]]:
    """
    For the given set of molecules produce a coverage report of the smirks parameters used in typing the molecule.
    Also try and produce charges for the molecules as some may have missing bccs.

    Parameters
    ----------
    input_molecules: The List of 3d sdf files to run the coverage on or the name of the directory containing the files.
    forcefield_name: The name of the openFF forcefield to run the coverage for.
    processors: The number of processors we can use to build the coverage report
    output_directory: The directory the files will be writen too if they pass the coverage step
    delete_existing: If any files currently exist remove them and run again

    Returns
    -------
        coverage_report: A dictionary split into parameter types which lists the number of occurrences of each parameter.
    """
    # make the file if needed
    try:
        os.makedirs(os.path.join(output_directory, "error_mols"))
    except OSError:
        if delete_existing:
            shutil.rmtree(output_directory)
            os.makedirs(os.path.join(output_directory, "error_mols"))
        else:
            raise Exception(f'Output directory {output_directory} already exists. '
                             'Specify `delete_existing=True` to remove.')

    if isinstance(input_molecules, list):
        input_files = input_molecules
    else:
        # check if its a dir then look for sdfs
        if os.path.isdir(input_molecules):
            input_files = glob.glob(os.path.join(input_molecules, "*.sdf"))
        else:
            input_files = [input_molecules, ]

    # now load each molecule they should already be unique
    molecules = [Molecule.from_file(mol_file, file_format="sdf", allow_undefined_stereo=True) for mol_file in input_files]

    coverage = {"Angles": {}, "Bonds": {}, "ProperTorsions": {}, "ImproperTorsions": {}, "vdW": {}}
    # a func to update the coverage dict
    def update_coverage(single_report):
        """
        Take the single molecule report and update the master coverage dict
        """
        for param_type in coverage.keys():
            for param_id, no in single_report[param_type].items():
                if param_id not in coverage[param_type]:
                    # set the param number if new
                    coverage[param_type][param_id] = no
                else:
                    # add onto the current parmeter count
                    coverage[param_type][param_id] += no

    # make the forcefield
    ff = ForceField(forcefield_name)
    # now run coverage on each molecule
    if processors is None or processors > 1:
        from multiprocessing import Pool

        with Pool(processes=processors) as pool:
            # generate the work list
            work_list = [pool.apply_async(single_molecule_coverage, (molecule, ff)) for molecule in molecules]
            for work in tqdm.tqdm(
                    work_list,
                    total=len(work_list),
                    ncols=80,
                    desc="{:30s}".format("Generating Coverage"),
            ):
                work = work.get()
                mol = work["molecule"]
                if "error" not in work:
                    update_coverage(work)
                    mol.to_file(os.path.join(output_directory, mol.name + ".sdf"), "sdf")
                else:
                    # Do we want the coverage for molecules that fail as well
                    mol.to_file(os.path.join(output_directory, "error_mols", mol.name + ".sdf"), "sdf")

    else:
        for molecule in tqdm.tqdm(molecules,
                                  total=len(molecules),
                                  ncols=80,
                                  desc="{:30s}".format("Generating Coverage")):
            report = single_molecule_coverage(molecule, ff)
            # update the master coverage dict
            mol: Molecule = report["molecule"]
            if "error" not in report:
                update_coverage(report)
                mol.to_file(os.path.join(output_directory, mol.name + ".sdf"), "sdf")
            else:
                # Do we want the coverage for molecules that fail as well
                mol.to_file(os.path.join(output_directory, "error_mols", mol.name + ".sdf"), "sdf")

    return coverage


def single_molecule_coverage(molecule: Molecule, forcefield: ForceField) -> Dict[str, Dict[str, int]]:
    """
    For a single molecule generate a coverage report and try to build an openmm system this will also highlight any missing parameters and dificulties with charging the molecule.

    Parameters
    ----------
    molecule: The openforcefield molecule object for which the report should be generated.
    ff: The openforcefield typing engine that should be used to check coverage and build an openmm system.

    Returns
    -------
        report: A dictionary of the coverage report along with any errors produced when trying to make an openmm system.
    """

    coverage = {"Angles": {}, "Bonds": {}, "ProperTorsions": {}, "ImproperTorsions": {}, "vdW": {}}
    labels = forcefield.label_molecules(molecule.to_topology())[0]
    for param_type, params in labels.items():
        for param in params.values():
            if param.id not in coverage[param_type]:
                coverage[param_type][param.id] = 1
            else:
                coverage[param_type][param.id] += 1
    # now generate a system this will catch any missing parameters
    # and molecules that can not be charged
    try:
        _ = forcefield.create_openmm_system(molecule.to_topology())
    except Exception as e:
        coverage["error"] = e

    coverage["molecule"] = molecule
    return coverage
