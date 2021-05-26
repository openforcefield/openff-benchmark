from openff.benchmark.utils.generate_conformers import generate_conformers
from openff.benchmark.cli import validate, cli
import inspect
import os
import pytest
import glob
import shutil
from openff.benchmark.utils.utils import get_data_file_path
from openff.toolkit.utils.toolkits import GLOBAL_TOOLKIT_REGISTRY, OpenEyeToolkitWrapper, RDKitToolkitWrapper
from openff.toolkit.topology import Molecule
from click.testing import CliRunner
runner = CliRunner()

def test_openeye_deregistered():
    for toolkitwrapper in GLOBAL_TOOLKIT_REGISTRY.registered_toolkits:
        assert not isinstance(toolkitwrapper, OpenEyeToolkitWrapper)
    to_smiles = GLOBAL_TOOLKIT_REGISTRY.resolve('to_smiles')
    assert not isinstance(to_smiles.__self__, OpenEyeToolkitWrapper)
    assert isinstance(to_smiles.__self__, RDKitToolkitWrapper)


def test_dont_overwrite_output_directory(tmpdir):
    with tmpdir.as_cwd():
        test_name = inspect.stack()[0].function
        input_dir = get_data_file_path('1-validate_and_assign_graphs_and_confs')
        output_dir = os.path.join(test_name, '2-generate_conformers')
        response = runner.invoke(cli, ["preprocess", "generate-conformers",
                                       "-o", output_dir,
                                       input_dir],
                                 catch_exceptions=False)
        with pytest.raises(Exception, match='Specify `--delete-existing` to remove'):
            response = runner.invoke(cli, ["preprocess", "generate-conformers",
                                           "-o", output_dir,
                                           input_dir],
                                     catch_exceptions=False)


# test loading a mix of graph and 3D molecules, generating up to 10 confs
def test_generate_conformers(tmpdir):
    with tmpdir.as_cwd():
        # test_name = inspect.stack()[0].function
        input_dir = get_data_file_path('1-validate_and_assign_graphs_and_confs')
        output_dir = '2-generate_conformers'

        # generate_conformers(input_dir, output_dir)
        response = runner.invoke(cli, ["preprocess", "generate-conformers",
                                       "-o", output_dir,
                                       input_dir],
                                 catch_exceptions=False)

        ## BBB-00000 starts with two conformers, so many more conformers should be created
        bbb0_confs = glob.glob(os.path.join(output_dir, 'BBB-00000-*.sdf'))
        assert len(bbb0_confs) > 3

        ## BBB-00001 starts with a one conformer, so many more conformers should be created
        bbb1_confs = glob.glob(os.path.join(output_dir, 'BBB-00001-*.sdf'))
        assert len(bbb1_confs) > 2

        ## BBB-00002 starts with one conformer.
        # It is rigid so only one conformer should be created
        bbb2_confs = glob.glob(os.path.join(output_dir, 'BBB-00002-*.sdf'))
        assert len(bbb2_confs) == 1

        ## BBB-00003 starts with 12 conformers.
        # We should see 12 output confs here, since we NEVER delete user confs
        bbb3_confs = glob.glob(os.path.join(output_dir, 'BBB-00003-*.sdf'))
        assert len(bbb3_confs) == 12

        ## BBB-00004 starts with a smi and NO conformers.
        # It is very flexible so 10 conformers should be output.
        #bbb4_confs = glob.glob(os.path.join(output_dir, 'BBB-00004-*.sdf'))
        #assert len(bbb4_confs) == 10


def test_generate_conformers_add(tmpdir):
    with tmpdir.as_cwd():
        # test_name = inspect.stack()[0].function
        input_dir = '1-validate_and_assign_graphs_and_confs'
        # Make a copy of this directory, otherwise we'll contaminate the original when we add a new mol for the test
        shutil.copytree(get_data_file_path(input_dir), input_dir)
        output_dir = '2-generate_conformers'
        # generate_conformers(input_dir, output_dir)
        response = runner.invoke(cli, ["preprocess", "generate-conformers",
                                       "-o", output_dir,
                                       input_dir],
                                 catch_exceptions=False)
        initial_confs = glob.glob(os.path.join(output_dir, 'BBB-*.sdf'))
        initial_confs = [os.path.basename(filename) for filename in initial_confs]

        # now add a new ridiculously flexible molecule to dir
        mol = Molecule.from_smiles('CCCCC[C@H](COCOC)COCCOCCCCCCC')
        mol.generate_conformers()
        mol.properties["group_name"] = "BBB"
        mol.properties["molecule_index"] = "99999"
        mol.to_file(os.path.join(input_dir, "BBB-99999-00.sdf"), "sdf")

        response = runner.invoke(cli, ["preprocess", "generate-conformers",
                                       "-o", output_dir,
                                       "--add",
                                       input_dir],
                                 catch_exceptions=False)

        final_confs = glob.glob(os.path.join(output_dir, 'BBB-*.sdf'))
        final_confs = [os.path.basename(filename) for filename in final_confs]
        assert 'BBB-99999-00.sdf' in final_confs
        assert 'BBB-99999-09.sdf' in final_confs
        assert len(final_confs) == len(initial_confs) + 10


def test_bad_macrocycle(tmpdir):
    with tmpdir.as_cwd():
        # test_name = inspect.stack()[0].function
        input_dir = get_data_file_path('1-validate_and_assign_graphs_and_confs_bad_macrocycle')
        output_dir = '2-generate_conformers'

        # generate_conformers(input_dir, output_dir)
        response = runner.invoke(cli, ["preprocess", "generate-conformers",
                                       "-o", output_dir,
                                       input_dir],
                                 catch_exceptions=False)

        # JAN_00203 has a macrocycle that RDKit generates bad conformers for. These conformers
        # have twisted double bonds and can't be parsed by subsequent processing steps.
        jan_203_confs = glob.glob(os.path.join(output_dir, 'JAN-00203-*.sdf'))
        assert len(jan_203_confs) == 1


def test_generate_conformers_add_bad_macrocycle(tmpdir):
    with tmpdir.as_cwd():
        # test_name = inspect.stack()[0].function
        input_dir = '1-validate_and_assign_graphs_and_confs'
        # Make a copy of this directory, otherwise we'll contaminate the original when we add a new mol for the test
        shutil.copytree(get_data_file_path(input_dir), input_dir)
        output_dir = '2-generate_conformers'
        # generate_conformers(input_dir, output_dir)
        response = runner.invoke(cli, ["preprocess", "generate-conformers",
                                       "-o", output_dir,
                                       input_dir],
                                 catch_exceptions=False)
        initial_confs = glob.glob(os.path.join(output_dir, 'BBB-*.sdf'))
        initial_confs = [os.path.basename(filename) for filename in initial_confs]

        # now add a new molecule to dir
        mol = Molecule.from_file(get_data_file_path("1-validate_and_assign_graphs_and_confs_bad_macrocycle/JAN-00203-00.sdf"))
        mol.properties["group_name"] = "BBB"
        mol.properties["molecule_index"] = "99999"
        mol.to_file(os.path.join(input_dir, "BBB-99999-00.sdf"), "sdf")

        response = runner.invoke(cli, ["preprocess", "generate-conformers",
                                       "-o", output_dir,
                                       "--add",
                                       input_dir],
                                 catch_exceptions=False)
        final_confs = glob.glob(os.path.join(output_dir, 'BBB-*.sdf'))
        final_confs = [os.path.basename(filename) for filename in final_confs]
        assert 'BBB-99999-00.sdf' in final_confs
        assert 'BBB-99999-01.sdf' not in final_confs
        assert len(final_confs) == len(initial_confs) + 1
        error_mols = glob.glob(os.path.join(output_dir, 'error_mols', '*'))
        assert len(error_mols) > 0
