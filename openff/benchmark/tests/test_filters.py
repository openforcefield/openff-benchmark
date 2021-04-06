from openff.benchmark.utils.filters import smirks_filter
from openff.benchmark.utils.utils import get_data_file_path
from openff.benchmark.cli import cli
from openforcefield.topology import Molecule
import glob
import os
import pytest
from click.testing import CliRunner
runner = CliRunner()


def test_single_smirks():
    """Test filtering based on a single smikrs."""
    molecules = []
    for i in [0, 1, 2, 3, 5]:
        molecules.append(Molecule.from_file(get_data_file_path(f'1-validate_and_assign_graphs_and_confs/BBB-0000{i}-00.sdf'), "sdf" ))
    # filter P should only be one molecule
    result = smirks_filter(input_molecules=molecules, filtered_smirks=["[P:1]"], processors=1)
    assert result.n_filtered == 1
    assert result.n_molecules == 4


def test_double_smirks():
    """Test filtering based on 2 different smirks patterns."""
    molecules = []
    for i in [0, 1, 2, 3, 5]:
        molecules.append(
            Molecule.from_file(get_data_file_path(f'1-validate_and_assign_graphs_and_confs/BBB-0000{i}-00.sdf'), "sdf"))
    # filter P should only be one molecule
    result = smirks_filter(input_molecules=molecules, filtered_smirks=["[P:1]", "[F:1]"], processors=1)
    assert result.n_filtered == 2
    assert result.n_molecules == 3


def test_cli_move_molecules(tmpdir):
    """Make sure that the cli can correctly move the molecules to the passed and fail directories."""
    with tmpdir.as_cwd():
        input_folder = get_data_file_path('1-validate_and_assign_graphs_and_confs')
        n_input_moles = len(glob.glob(os.path.join(input_folder, "*.sdf")))
        test_dir = '5-smirks_filter'
        response = runner.invoke(cli, ["filter", "smirks",
                                       input_folder,
                                       test_dir,
                                       "-p", 1,
                                       "-s", "[P:1]"],
                                 catch_exceptions=False)
        n_out_mols = len(glob.glob(os.path.join(test_dir, "*.sdf")))
        # this should only remove 1 molecule with 2 conformers
        assert n_out_mols == n_input_moles - 2
        n_error_mols = len(glob.glob(os.path.join(test_dir, "error_mols", "*.sdf")))
        assert n_error_mols == 2