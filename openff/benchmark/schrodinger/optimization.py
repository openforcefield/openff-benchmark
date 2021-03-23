"""
Components for energy minimization / geometry optimization of molecules for FF benchmarking.

These components leverage the QCArchive software stack for execution.

"""

import os
import shutil
import logging
import subprocess
import warnings
from collections import defaultdict

# DO NOT USE; MESSES UP LOGGING
# consider overriding the loglevel of any particularly noisy loggers individual
# logging.disable(logging.WARNING)

from openforcefield.topology import Molecule
from ..utils.io import mols_from_paths


def sdf_to_mae(sdf_file, mae_file, schrodinger_path, append=True):
    command = [
        os.path.join(schrodinger_path, "utilities", "sdconvert"),
        "-isd",
        sdf_file,
        "-omae",
        mae_file,
    ]
    if append:
        command.append("-a")
    subprocess.run(command)


def mae_to_sdf(mae_file, sdf_file, schrodinger_path, append=True):
    command = [
        os.path.join(schrodinger_path, "utilities", "sdconvert"),
        "-imae",
        mae_file,
        "-osd",
        sdf_file,
    ]
    if append:
        command.append("-a")
    subprocess.run(command)


def optimization(
    input_paths,
    schrodinger_path,
    host_settings,
    opls_dir,
    recursive,
    output_path,
    delete_existing=False,
):
    """Submit SDF molecules from given directory to a Schrodinger macromodel OPLS3e minimization

    Parameters
    ----------
    input_paths : str or list of str
        Input SDF molecules or paths containing SDFs
    schrodinger_path : str
        Path to Schrodinger executables
    host_settings : str
        Settings for Schrodinger job of format '<host>:<max number of jobs>'
    opls_dir : str
        Path to custom OPLS3e parameter set
    recursive : bool
        If True, recursively load SDFs from any directories given in `input_paths`.
    output_path: str
        Path where output sdf files should be written to
    delete_existing : bool (False)
            If True, delete existing directory if present.
    """
    try:
        os.makedirs(output_path)
    except OSError:
        if delete_existing:
            shutil.rmtree(output_path)
            os.makedirs(output_path)
        else:
            raise Exception(
                f"Output directory {output_path} already exists. "
                "Specify `delete_existing=True` to remove or select a different output directory."
            )

    mae_file = os.path.join(output_path, "mmod_input.maegz")
    sdf_file = os.path.join(output_path, "mmod_input.sdf")

    mols = mols_from_paths(input_paths, recursive=recursive)
    if len(mols) == 0:
        raise Exception("No molecules are selected for the ffbuilder job. "
                        "Check the input path.")
        
    with open(sdf_file, "w") as file:
        for i, mol in enumerate(mols):
            if opls_dir is None:
                mol.properties["method"] = "opls3e_default"
                mol.properties["basis"] = "opls"
                mol.properties["program"] = "macromodel"
            else:
                mol.properties["method"] = "opls3e_custom"
                mol.properties["basis"] = "opls"
                mol.properties["program"] = "macromodel"
            mol.to_file(file, "SDF")

    sdf_to_mae(sdf_file, mae_file, schrodinger_path)

    com_file = os.path.join(output_path, "mmod_openff_benchmark.com")
    with open(com_file, "w") as outfile:
        outfile.write("mmod_input.maegz\n")
        outfile.write("mmod_output.maegz\n")
        outfile.write(
            " MMOD       0      1      0      0     0.0000     0.0000     0.0000     0.0000\n"
            " DEBG      55      0      0      0     0.0000     0.0000     0.0000     0.0000\n"
            " FFLD      16      1      0      0     1.0000     0.0000     0.0000     0.0000\n"
            " BDCO       0      0      0      0    41.5692 99999.0000     0.0000     0.0000\n"
            " CRMS       0      0      0      0     0.0000     0.5000     0.0000     0.0000\n"
            " BGIN       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n"
            " READ       0      0      0      0     0.0000     0.0000     0.0000     0.0000\n"
            " CONV       1      0      0      0     5.0e-9     0.0000     0.0000     0.0000\n"
            " MINI      10      0   1500      0     0.0000     0.0000     0.0000     0.0000\n"
            " END        0      0      0      0     0.0000     0.0000     0.0000     0.0000\n"
        )

    command = [
        os.path.join(schrodinger_path, "macromodel"),
        "mmod_openff_benchmark.com",
    ]
    if host_settings is not None:
        command.append("-HOST")
        command.append(f"{host_settings}")
    if opls_dir is not None:
        command.append("-OPLSDIR")
        command.append(os.path.abspath(opls_dir))
    subprocess.run(command, cwd=os.path.join(output_path))


def postprocess(
    input_paths,
    schrodinger_path,
    output_path,
    delete_existing=False,
    keep_existing=True,
):
    """Postprocess output from Schrodinger macromodel OPLS3e minimization.

    Parameters
    ----------
    input_paths : str or list of str
        Schrodinger *.mae or *.maegz file (output from schrodinger optimize).
    schrodinger_path : str
        Path to Schrodinger executables
    output_path: str
        Path where output sdf files should be written to
    delete_existing : bool (False)
            If True, delete existing directory if present.
    keep_existing : bool (True)
            If True, keep existing files in export directory.
            Relies *only* on filepaths of existing files for determining match.
    """

    mols = []
    for path in input_paths:
        if os.path.isfile(path) and (
            path.split(".")[-1].lower() == "maegz"
            or path.split(".")[-1].lower() == "mae"
        ):
            sdf_file = ".".join(path.split(".")[:-1] + ["sdf"])
            mae_to_sdf(path, sdf_file, schrodinger_path, append=False)
            tmp_mols = Molecule.from_file(sdf_file, "SDF")
            if type(tmp_mols) == Molecule:
                mols.append(tmp_mols)
            else:
                mols += tmp_mols

    methods = set()
    for mol in mols:
        mol.properties.pop("initial_energy")
        mol.properties[
            "final_energy"
        ] = f'{mol.properties["r_mmod_Potential_Energy-OPLS3e"]} kJ / mole'
        mol.properties["molecule_index"] = f"{mol.properties['molecule_index']:05d}"
        mol.properties["conformer_index"] = f"{mol.properties['conformer_index']:02d}"
        methods.add(mol.properties["method"])

    os.makedirs(output_path, exist_ok=True)
    for method in methods:
        output_directory = os.path.join(output_path, method)
        try:
            os.makedirs(output_directory)
        except OSError:
            if delete_existing:
                shutil.rmtree(output_directory)
                os.makedirs(output_directory)
            elif keep_existing:
                warnings.warn(
                    f"Output directory {output_directory} exists and will be used for export; existing data files will not be replaced."
                    "Specify `delete_existing=True` to replace."
                )
                pass
            else:
                raise Exception(
                    f"Output directory {output_directory} already exists. "
                    "Specify `delete_existing=True` to remove or specify another output directory."
                )
    for mol in mols:
        outfile = os.path.join(output_path, mol.properties["method"], f"{mol.name}.sdf")
        if (not delete_existing) and os.path.exists(outfile):
            print(f"... '{mol.name}' : skipping SDF exists")
            continue
        mol.to_file(outfile, "SDF")
