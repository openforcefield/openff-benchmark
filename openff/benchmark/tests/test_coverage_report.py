from openff.benchmark.utils.coverage_report import generate_coverage_report
from openff.benchmark.utils.utils import get_data_file_path
from openff.benchmark.cli import cli
from openforcefield.topology import Molecule
import json
import glob
import os
import pytest
from click.testing import CliRunner
runner = CliRunner()

# Test single file
@pytest.mark.parametrize('n_procs', [1,2])
def test_single_file(n_procs):
    molecules = Molecule.from_file(get_data_file_path('1-validate_and_assign_graphs_and_confs/BBB-00000-00.sdf'), "sdf" )
    coverage, success_mols, error_mols = generate_coverage_report(input_molecules=molecules,
                                                                  forcefield_name='openff_unconstrained-1.3.0.offxml',
                                                                  processors=n_procs)
    assert len(success_mols) == 1, error_mols[0][1]
    assert len(error_mols) == 0

# Test multiple files
@pytest.mark.parametrize('n_procs', [1,2])
def test_multiple_files(n_procs):
    molecules = [Molecule.from_file(get_data_file_path('1-validate_and_assign_graphs_and_confs/BBB-00000-00.sdf'), "sdf"),
               Molecule.from_file(get_data_file_path('1-validate_and_assign_graphs_and_confs/BBB-00001-00.sdf'), "sdf")]
    coverage, success_mols, error_mols = generate_coverage_report(input_molecules=molecules,
                                                                  forcefield_name='openff_unconstrained-1.3.0.offxml',
                                                                  processors=n_procs)
    assert len(success_mols) == 2
    assert len(error_mols) == 0
    assert coverage["passed_unique_molecules"] == 2
    assert coverage["total_unique_molecules"] == 2
    #raise Exception(error_mols)



# Test catching uncovered params
@pytest.mark.parametrize('n_procs',[1,2])
def test_error_missing_valence_param(n_procs):
    molecules = Molecule.from_file(get_data_file_path('missing_valence_params.sdf'), "sdf")
    coverage, success_mols, error_mols = generate_coverage_report(input_molecules=molecules,
                                                                  forcefield_name='openff_unconstrained-1.3.0.offxml',
                                                                  processors=n_procs)
    assert len(success_mols) == 0
    assert len(error_mols) == 1
    assert coverage["passed_unique_molecules"] == 0
    assert coverage["total_unique_molecules"] == 1
    assert "BondHandler was not able to find parameters" in str(error_mols[0][1])

# Test for a molecule covered by OpenFF but that causes antechamber to fail
# I can't actually come up with such a molecule -- Suggestions welcome
# Also note that charge gen may be getting skipped --
# Search "deregister_parameter_handler" in coverage_report.py
@pytest.mark.skip
def test_error_uncovered_antechamber_param():
    molecule = Molecule.from_file(get_data_file_path('sodium_carbide.sdf'), "sdf")
    coverage, success_mols, error_mols = generate_coverage_report(input_molecules=molecule,
                                                                  forcefield_name='openff_unconstrained-1.3.0.offxml')
    assert len(success_mols) == 0
    assert len(error_mols) == 1
    assert coverage["passed_unique_molecules"] == 0
    assert coverage["total_unique_molecules"] == 1
    assert "Command '['antechamber'" in str(error_mols[0][1])


def test_cli_move_all_confs(tmpdir):
    """
    Make sure that if a molecule passes all conformers are also moved.
    """

    with tmpdir.as_cwd():
        test_dir = '3-coverage_report'
        input_folder = get_data_file_path('1-validate_and_assign_graphs_and_confs')
        # get the number of input molecules and conformers
        n_input_moles = len(glob.glob(os.path.join(input_folder, "*.sdf")))
        response = runner.invoke(cli, ["preprocess", "coverage-report",
                                       "-p", 1,
                                       "-o", test_dir,
                                       input_folder],
                                 catch_exceptions=False)
        n_out_mols = len(glob.glob(os.path.join(test_dir, "*.sdf")))
        # assuming no molecules fail
        assert n_input_moles == n_out_mols
        n_error_mols = len(glob.glob(os.path.join(test_dir, "error_mols", "*.sdf")))
        assert n_error_mols == 0
        # load the coverage report and make sure the unique mols is correct
        with open(os.path.join(test_dir, "coverage_report.json"), "r") as data:
            report = json.load(data)

        assert report["passed_unique_molecules"] == 5
        assert report["total_unique_molecules"] == 5


def test_cli_error_mol(tmpdir):
    """Make sure that molecules that fail are correctly put in error mols folder."""

    with tmpdir.as_cwd():
        test_dir = '3-coverage_report'
        input_folder = "2-generate_conformers"
        os.mkdir(input_folder)
        mol = Molecule.from_file(get_data_file_path('missing_valence_params.sdf'), "sdf", allow_undefined_stereo=True)
        mol.properties["group_name"] = "OFF"
        mol.properties["molecule_index"] = "00000"
        mol.to_file(os.path.join(input_folder, "OFF-00000-00.sdf"), "sdf")
        response = runner.invoke(cli, ["preprocess", "coverage-report",
                                       "-p", 1,
                                       "-o", test_dir,
                                       input_folder],
                                 catch_exceptions=False)

        n_out_mols = len(glob.glob(os.path.join(test_dir, "*.sdf")))
        assert n_out_mols == 0
        n_error_mols = len(glob.glob(os.path.join(test_dir, "error_mols", "*.sdf")))
        assert n_error_mols == 1
        # load the coverage report and make sure the unique mols is correct
        with open(os.path.join(test_dir, "coverage_report.json"), "r") as data:
            report = json.load(data)

        assert report["passed_unique_molecules"] == 0
        assert report["total_unique_molecules"] == 1