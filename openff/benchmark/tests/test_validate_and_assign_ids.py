from openforcefield.utils.toolkits import GLOBAL_TOOLKIT_REGISTRY, OpenEyeToolkitWrapper, RDKitToolkitWrapper
from openff.benchmark.utils.validate_and_assign_ids import validate_and_assign
from openff.benchmark.cli import validate, cli
from openff.benchmark.utils.utils import get_data_file_path
import pytest
import os
import shutil
import inspect
import glob
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
        test_dir = '1-validate_and_assign_graphs_and_confs'
        input_mols = [get_data_file_path('input_single_mol_rigid.sdf')]
        response = runner.invoke(cli, ["preprocess", "validate",
                                       "-g", "BBB",
                                       "-o", test_dir,
                                       *input_mols],
                                 catch_exceptions=False)

        with pytest.raises(Exception, match='Specify `delete_existing=True` to remove'):
            response = runner.invoke(cli, ["preprocess", "validate",
                                           "-g", "BBB",
                                           "-o", test_dir,
                                           *input_mols],
                                     catch_exceptions=False)

def test_do_overwrite_output_directory(tmpdir):
    with tmpdir.as_cwd():
        test_dir = '1-validate_and_assign_graphs_and_confs'
        input_mols = [get_data_file_path('input_single_mol_rigid.sdf')]
        response = runner.invoke(cli, ["preprocess", "validate",
                                       "-g", "BBB",
                                       "-o", test_dir,
                                       *input_mols],
                                 catch_exceptions=False)

        response = runner.invoke(cli, ["preprocess", "validate",
                                       "-g", "BBB",
                                       "-o", test_dir,
                                       "--delete-existing",
                                       *input_mols],
                                 catch_exceptions=False)
class TestCLI:
    # single file single mol
    def test_single_file_single_mol(self, tmpdir):
        with tmpdir.as_cwd():
            test_dir = os.path.join('1-validate_and_assign_graphs_and_confs')

            input_mols = [get_data_file_path('input_single_mol_rigid.sdf')]
            response = runner.invoke(cli, ["preprocess", "validate",
                                           "-g", "BBB",
                                           "-o", test_dir,
                                           *input_mols],
                                     catch_exceptions=False)

            output_files = glob.glob(os.path.join(test_dir, '*.sdf'))
            output_files = [os.path.basename(fname) for fname in output_files]
            assert 'BBB-00000-00.sdf' in output_files
            assert len(output_files) == 1
            file_text = open(os.path.join(test_dir, 'BBB-00000-00.sdf')).read()
            assert """
>  <group_name>  (1) 
BBB""" in file_text
            assert """
>  <molecule_index>  (1) 
0""" in file_text
            assert """
>  <conformer_index>  (1) 
0""" in file_text

    # single file multi mol
    def test_single_file_multi_mol(self, tmpdir):
        with tmpdir.as_cwd():
            test_dir = '1-validate_and_assign_graphs_and_confs'
            if os.path.exists(test_dir):
                shutil.rmtree(test_dir)
            input_mols = [get_data_file_path('input_one_stereoisomer_and_multi_conf_flexible.sdf')]
            response = runner.invoke(cli, ["preprocess", "validate",
                                           "-g", "BBB",
                                           "-o", test_dir,
                                           *input_mols],
                                     catch_exceptions=False)

            output_files = glob.glob(os.path.join(test_dir, '*.sdf'))
            output_files = [os.path.basename(fname) for fname in output_files]

            assert sorted(output_files) == ['BBB-00000-00.sdf',
                                            'BBB-00001-00.sdf',
                                            'BBB-00001-01.sdf',
                                            'BBB-00001-02.sdf',
                                            'BBB-00001-03.sdf',
                                            'BBB-00001-04.sdf']

    # single file multi conformer
    @pytest.mark.xfail
    def test_single_file_multi_conformer(self, tmpdir):
        with tmpdir.as_cwd():
            test_dir = '1-validate_and_assign_graphs_and_confs'
            if os.path.exists(test_dir):
                shutil.rmtree(test_dir)
            input_mols = [get_data_file_path('input_multi_conf_flexible.sdf')]
            response = runner.invoke(cli, ["preprocess", "validate",
                                           "-g", "BBB",
                                           "-o", test_dir,
                                           *input_mols],
                                     catch_exceptions=False)
            output_files = glob.glob(os.path.join(test_dir, '*.sdf'))
            output_files = [os.path.basename(fname) for fname in output_files]
            assert 'BBB-00000-00.sdf' in output_files
            assert 'BBB-00000-09.sdf' in output_files
            assert len(output_files) == 10
        
    # multi file single mol
    def test_multi_file_single_mol(self, tmpdir):
        with tmpdir.as_cwd():
            test_dir =  '1-validate_and_assign_graphs_and_confs'
            #test_dir = 'test_multi_molecule_multi_sdf_graph_input'
            if os.path.exists(test_dir):
                shutil.rmtree(test_dir)
            #os.makedirs(test_dir)
            input_mols = [get_data_file_path('input_single_mol_flexible.sdf'),
                          get_data_file_path('input_one_stereoisomer.sdf')]
            response = runner.invoke(cli, ["preprocess", "validate",
                                           "-g", "BBB",
                                           "-o", test_dir,
                                           *input_mols],
                                     catch_exceptions=False)
            output_files = glob.glob(os.path.join(test_dir, '*.sdf'))
            output_files = [os.path.basename(fname) for fname in output_files]
            #assert 'BBB-00000.smi' in output_files
            #assert 'BBB-00001.smi' in output_files
            assert 'BBB-00000-00.sdf' in output_files
            assert 'BBB-00001-00.sdf' in output_files
            assert len(output_files) == 2
        
    # multi file multi mol
    def test_multi_file_multi_mol(self, tmpdir):
        with tmpdir.as_cwd():
            test_dir = '1-validate_and_assign_graphs_and_confs'
            #test_dir = 'test_multi_molecule_multi_sdf_graph_input'
            if os.path.exists(test_dir):
                shutil.rmtree(test_dir)
            #os.makedirs(test_dir)
            input_mols = [get_data_file_path('input_multi_conf_flexible.sdf'),
                          get_data_file_path('input_all_stereoisomers.sdf')]
            response = runner.invoke(cli, ["preprocess", "validate",
                                           "-g", "BBB",
                                           "-o", test_dir,
                                           *input_mols],
                                     catch_exceptions=False)
            output_files = glob.glob(os.path.join(test_dir, '*.sdf'))
            output_files = [os.path.basename(fname) for fname in output_files]
            #assert 'BBB-00000.smi' in output_files
            #assert 'BBB-00011.smi' in output_files
            #assert 'BBB-00012.smi' not in output_files
            assert 'BBB-00000-00.sdf' in output_files
            assert 'BBB-00011-00.sdf' in output_files
            assert len(output_files) == 12

    # multi file multi mol duplicates (ok)
    def test_multi_file_multi_mol_duplicates(self, tmpdir):
        with tmpdir.as_cwd():
            # raise Exception(inspect.stack()[0])
            test_dir = '1-validate_and_assign_graphs_and_confs'
            # test_dir = 'test_multi_molecule_multi_sdf_graph_input'
            if os.path.exists(test_dir):
                shutil.rmtree(test_dir)
            # os.makedirs(test_dir)
            input_mols = [get_data_file_path('input_one_stereoisomer.sdf'),
                          get_data_file_path('input_all_stereoisomers.sdf')]
            response = runner.invoke(cli, ["preprocess", "validate",
                                           "-g", "BBB",
                                           "-o", test_dir,
                                           *input_mols],
                                     catch_exceptions=False)
            output_files = glob.glob(os.path.join(test_dir, '*.sdf'))
            output_files = [os.path.basename(fname) for fname in output_files]
            #assert 'BBB-00000.smi' in output_files
            #assert 'BBB-00007.smi' in output_files
            assert 'BBB-00000-00.sdf' in output_files
            assert 'BBB-00000-01.sdf' not in output_files
            assert 'BBB-00007-00.sdf' in output_files
            assert 'BBB-00008-00.sdf' not in output_files
            assert len(output_files) == 8

    def test_add(self, tmpdir):

        with tmpdir.as_cwd():
            #raise Exception(inspect.stack()[0])
            test_dir = '1-validate_and_assign_graphs_and_confs'
            #test_dir = 'test_multi_molecule_multi_sdf_graph_input'
            if os.path.exists(test_dir):
                shutil.rmtree(test_dir)
            #os.makedirs(test_dir)
            input_mols = [get_data_file_path('input_one_stereoisomer.sdf')]
            response = runner.invoke(cli, ["preprocess", "validate",
                                           "-g", "BBB",
                                           "-o", test_dir,
                                           *input_mols],
                                     catch_exceptions=False)
            input_mols2 = [get_data_file_path('input_all_stereoisomers.sdf')]
            response = runner.invoke(cli, ["preprocess", "validate",
                                           "--add",
                                           "-g", "BBB",
                                           "-o", test_dir,
                                           *input_mols2],
                                     catch_exceptions=False)

            output_files = glob.glob(os.path.join(test_dir, '*.sdf'))
            output_files = [os.path.basename(fname) for fname in output_files]
            error_files = glob.glob(os.path.join(test_dir, 'error_mols', '*.sdf'))
            error_files = [os.path.basename(fname) for fname in error_files]
