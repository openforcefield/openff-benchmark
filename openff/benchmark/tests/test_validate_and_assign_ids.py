from openforcefield.utils.toolkits import GLOBAL_TOOLKIT_REGISTRY, OpenEyeToolkitWrapper, RDKitToolkitWrapper
from openff.benchmark.utils.validate_and_assign_ids import validate_and_assign
from openff.benchmark.cli import validate, cli
from openff.benchmark.utils.utils import get_data_file_path
import pytest
import os
import shutil
import csv
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
        test_dir = '1-validate_and_assign'
        input_mols = [get_data_file_path('input_single_mol_rigid.sdf')]
        response = runner.invoke(cli, ["preprocess", "validate",
                                       "-g", "BBB",
                                       "-o", test_dir,
                                       *input_mols],
                                 catch_exceptions=False)

        with pytest.raises(Exception, match='Specify `--delete-existing` to remove'):
            response = runner.invoke(cli, ["preprocess", "validate",
                                           "-g", "BBB",
                                           "-o", test_dir,
                                           *input_mols],
                                     catch_exceptions=False)

def test_do_overwrite_output_directory(tmpdir):
    with tmpdir.as_cwd():
        test_dir = '1-validate_and_assign'
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


def test_add_and_delete_existing_error(tmpdir):
    with tmpdir.as_cwd():
        test_dir = '1-validate_and_assign'
        input_mols = [get_data_file_path('input_one_stereoisomer.sdf')]
        with pytest.raises(Exception, match='Can not specify BOTH') as context:
            response = runner.invoke(cli, ["preprocess", "validate",
                                           "--add",
                                           "--delete-existing",
                                           "-g", "BBB",
                                           "-o", test_dir,
                                           *input_mols],
                                     catch_exceptions=False)

class TestCLI:
    # single file single mol
    def test_single_file_single_mol(self, tmpdir):
        with tmpdir.as_cwd():
            test_dir = '1-validate_and_assign'

            input_mols = [get_data_file_path('input_single_mol_rigid.sdf')]
            input_mols = [os.path.abspath(input_mol) for input_mol in input_mols]
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

            error_files = glob.glob(os.path.join(test_dir, 'error_mols', '*.sdf'))
            error_files = [os.path.basename(fname) for fname in error_files]
            assert error_files == []



    # single file multi mol
    def test_single_file_multi_mol(self, tmpdir):
        with tmpdir.as_cwd():
            test_dir = '1-validate_and_assign'
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

            error_files = glob.glob(os.path.join(test_dir, 'error_mols', '*.sdf'))
            error_files = [os.path.basename(fname) for fname in error_files]
            assert error_files == []



    # single file multi conformer
    def test_single_file_multi_conformer(self, tmpdir):
        with tmpdir.as_cwd():
            test_dir = '1-validate_and_assign'
            input_mols = [get_data_file_path('input_five_confs_flexible.sdf')]
            response = runner.invoke(cli, ["preprocess", "validate",
                                           "-g", "BBB",
                                           "-o", test_dir,
                                           *input_mols],
                                     catch_exceptions=False)
            output_files = glob.glob(os.path.join(test_dir, '*.sdf'))
            output_files = [os.path.basename(fname) for fname in output_files]
            assert sorted(output_files) == ['BBB-00000-00.sdf',
                                            'BBB-00000-01.sdf',
                                            'BBB-00000-02.sdf',
                                            'BBB-00000-03.sdf',
                                            'BBB-00000-04.sdf',
                                            ]
            error_files = glob.glob(os.path.join(test_dir, 'error_mols', '*.sdf'))
            error_files = [os.path.basename(fname) for fname in error_files]
            assert error_files == []


    # multi file single mol
    def test_multi_file_single_mol(self, tmpdir):
        with tmpdir.as_cwd():
            test_dir = '1-validate_and_assign'
            input_mols = [get_data_file_path('input_single_mol_flexible.sdf'),
                          get_data_file_path('input_one_stereoisomer.sdf')]
            response = runner.invoke(cli, ["preprocess", "validate",
                                           "-g", "BBB",
                                           "-o", test_dir,
                                           *input_mols],
                                     catch_exceptions=False)
            output_files = glob.glob(os.path.join(test_dir, '*.sdf'))
            output_files = [os.path.basename(fname) for fname in output_files]
            assert 'BBB-00000-00.sdf' in output_files
            assert 'BBB-00001-00.sdf' in output_files
            assert len(output_files) == 2
            error_files = glob.glob(os.path.join(test_dir, 'error_mols', '*.sdf'))
            error_files = [os.path.basename(fname) for fname in error_files]
            assert error_files == []


    # multi file single mol
    def test_multi_file_single_mol_redundant_conf(self, tmpdir):
        with tmpdir.as_cwd():
            test_dir =  '1-validate_and_assign'
            input_mols = [get_data_file_path('input_single_mol_rigid.sdf'),
                          get_data_file_path('input_single_mol_rigid_translated.sdf')]
            response = runner.invoke(cli, ["preprocess", "validate",
                                           "-g", "BBB",
                                           "-o", test_dir,
                                           *input_mols],
                                     catch_exceptions=False)
            output_files = glob.glob(os.path.join(test_dir, '*.sdf'))
            output_files = [os.path.basename(fname) for fname in output_files]
            assert sorted(output_files) == ['BBB-00000-00.sdf']
            error_files = glob.glob(os.path.join(test_dir, 'error_mols', '*.sdf'))
            error_files = [os.path.basename(fname) for fname in error_files]
            assert error_files == ['error_mol_0.sdf']
            error_txts = glob.glob(os.path.join(test_dir, 'error_mols', '*.txt'))
            error_txts = [open(fname).read() for fname in error_txts]
            assert 'Duplicate molecule conformer input detected' in error_txts[0]


    # multi file multi mol
    def test_multi_file_multi_mol(self, tmpdir):
        with tmpdir.as_cwd():
            test_dir = '1-validate_and_assign'
            input_mols = [get_data_file_path('input_five_confs_flexible.sdf'),
                          get_data_file_path('input_eight_stereoisomers.sdf')
                          ]
            response = runner.invoke(cli, ["preprocess", "validate",
                                           "-g", "BBB",
                                           "-o", test_dir,
                                           *input_mols],
                                     catch_exceptions=False)
            output_files = glob.glob(os.path.join(test_dir, '*.sdf'))
            output_files = [os.path.basename(fname) for fname in output_files]
            assert sorted(output_files) == ['BBB-00000-00.sdf', # The first input has 5 confs of the same molecule
                                            'BBB-00000-01.sdf',
                                            'BBB-00000-02.sdf',
                                            'BBB-00000-03.sdf',
                                            'BBB-00000-04.sdf',
                                            'BBB-00001-00.sdf',# The there are 8 different stereoisomers with 1 conf each
                                            'BBB-00002-00.sdf',
                                            'BBB-00003-00.sdf',
                                            'BBB-00004-00.sdf',
                                            'BBB-00005-00.sdf',
                                            'BBB-00006-00.sdf',
                                            'BBB-00007-00.sdf',
                                            'BBB-00008-00.sdf',
                                            ]

            error_files = glob.glob(os.path.join(test_dir, 'error_mols', '*.sdf'))
            error_files = [os.path.basename(fname) for fname in error_files]
            assert error_files == []

    # multi file multi mol duplicates (ok)
    def test_multi_file_multi_mol_duplicates(self, tmpdir):
        with tmpdir.as_cwd():
            test_dir = '1-validate_and_assign'
            input_mols = [get_data_file_path('input_one_stereoisomer.sdf'),
                          get_data_file_path('input_eight_stereoisomers.sdf')]
            response = runner.invoke(cli, ["preprocess", "validate",
                                           "-g", "BBB",
                                           "-o", test_dir,
                                           *input_mols],
                                     catch_exceptions=False)
            output_files = glob.glob(os.path.join(test_dir, '*.sdf'))
            output_files = [os.path.basename(fname) for fname in output_files]
            assert sorted(output_files) == ['BBB-00000-00.sdf',
                                            'BBB-00001-00.sdf',
                                            'BBB-00002-00.sdf',
                                            'BBB-00003-00.sdf',
                                            'BBB-00004-00.sdf',
                                            'BBB-00005-00.sdf',
                                            'BBB-00006-00.sdf',
                                            'BBB-00007-00.sdf',
                                            ]
            error_files = glob.glob(os.path.join(test_dir, 'error_mols', '*.sdf'))
            assert len(error_files) == 1

    def test_add(self, tmpdir):

        with tmpdir.as_cwd():
            test_dir = '1-validate_and_assign'
            input_mols = [get_data_file_path('input_one_stereoisomer.sdf')]
            response = runner.invoke(cli, ["preprocess", "validate",
                                           "-g", "BBB",
                                           "-o", test_dir,
                                           *input_mols],
                                     catch_exceptions=False)

            output_files = glob.glob(os.path.join(test_dir, '*.sdf'))
            output_files = [os.path.basename(fname) for fname in output_files]
            assert sorted(output_files) == ['BBB-00000-00.sdf']

            error_files = glob.glob(os.path.join(test_dir, 'error_mols', '*.sdf'))
            error_files = [os.path.basename(fname) for fname in error_files]
            assert error_files == []

            # Test that output names were correctly assigned
            output_name_assignments = []
            with open(os.path.join(test_dir, 'name_assignments.csv')) as of:
                csv_reader = csv.reader(of)
                for row in csv_reader:
                    new_row = [os.path.basename(i) for i in row]
                    output_name_assignments.append(new_row)
            assert output_name_assignments == [['orig_name', 'orig_file', 'orig_file_index', 'out_file_name'],
                                               ['PDB_DB00136_99', 'input_one_stereoisomer.sdf', '0', 'BBB-00000-00'],
                                              ]

            # Run again, with a partially overlapping set of molecules
            input_mols2 = [get_data_file_path('input_eight_stereoisomers.sdf')]
            response = runner.invoke(cli, ["preprocess", "validate",
                                           "--add",
                                           "-g", "BBB",
                                           "-o", test_dir,
                                           *input_mols2],
                                     catch_exceptions=False)

            output_files = glob.glob(os.path.join(test_dir, '*.sdf'))
            output_files = [os.path.basename(fname) for fname in output_files]
            assert sorted(output_files) == ['BBB-00000-00.sdf',
                                            'BBB-00001-00.sdf',
                                            'BBB-00002-00.sdf',
                                            'BBB-00003-00.sdf',
                                            'BBB-00004-00.sdf',
                                            'BBB-00005-00.sdf',
                                            'BBB-00006-00.sdf',
                                            'BBB-00007-00.sdf',
                                            ]

            error_files = glob.glob(os.path.join(test_dir, 'error_mols', '*.sdf'))
            error_files = [os.path.basename(fname) for fname in error_files]
            assert error_files == ['error_mol_0.sdf']
            error_txts = glob.glob(os.path.join(test_dir, 'error_mols', '*.txt'))
            error_txts = [open(fname).read() for fname in error_txts]
            assert "Input molecule graph is already present in output" in error_txts[0]

            # Test that output names were correctly assigned
            output_name_assignments = []
            with open(os.path.join(test_dir, 'name_assignments.csv')) as of:
                csv_reader = csv.reader(of)
                for row in csv_reader:
                    new_row = [os.path.basename(i) for i in row]
                    output_name_assignments.append(new_row)
            assert output_name_assignments == [['orig_name', 'orig_file', 'orig_file_index', 'out_file_name'],
                                               ['PDB_DB00136_99', 'input_one_stereoisomer.sdf', '0', 'BBB-00000-00'],
                                               ['PDB_DB00136_01', 'input_eight_stereoisomers.sdf', '1', 'BBB-00001-00'],
                                               ['PDB_DB00136_02', 'input_eight_stereoisomers.sdf', '2', 'BBB-00002-00'],
                                               ['PDB_DB00136_03', 'input_eight_stereoisomers.sdf', '3', 'BBB-00003-00'],
                                               ['PDB_DB00136_04', 'input_eight_stereoisomers.sdf', '4', 'BBB-00004-00'],
                                               ['PDB_DB00136_05', 'input_eight_stereoisomers.sdf', '5', 'BBB-00005-00'],
                                               ['PDB_DB00136_06', 'input_eight_stereoisomers.sdf', '6', 'BBB-00006-00'],
                                               ['PDB_DB00136_07', 'input_eight_stereoisomers.sdf', '7', 'BBB-00007-00'],
                                              ]


    def test_add_doesnt_overwrite_error_mols(self, tmpdir):
        """
        Run add multiple times, such that error mols are generated by two separate invocations. Then, make sure
        that the error outputs don't overwrite each other
        """
        with tmpdir.as_cwd():
            test_dir = '1-validate_and_assign'
            input_mols = [get_data_file_path('input_single_mol_rigid.sdf')]
            response = runner.invoke(cli, ["preprocess", "validate",
                                           "-g", "BBB",
                                           "-o", test_dir,
                                           *input_mols],
                                     catch_exceptions=False)

            response = runner.invoke(cli, ["preprocess", "validate",
                                           '--add',
                                           "-g", "BBB",
                                           "-o", test_dir,
                                           *input_mols],
                                     catch_exceptions=False)

            response = runner.invoke(cli, ["preprocess", "validate",
                                           '--add',
                                           "-g", "BBB",
                                           "-o", test_dir,
                                           *input_mols],
                                     catch_exceptions=False)
            output_files = glob.glob(os.path.join(test_dir, '*.sdf'))
            output_files = [os.path.basename(fname) for fname in output_files]
            assert sorted(output_files) == ['BBB-00000-00.sdf']
            error_files = glob.glob(os.path.join(test_dir, 'error_mols', '*.sdf'))
            error_files = [os.path.basename(fname) for fname in error_files]
            assert sorted(error_files )== ['error_mol_0.sdf', 'error_mol_1.sdf']
            error_txts = sorted(glob.glob(os.path.join(test_dir, 'error_mols', '*.txt')))
            error_txts = [open(fname).read() for fname in error_txts]
            assert "Input molecule graph is already present in output" in error_txts[0]
            assert "Input molecule graph is already present in output" in error_txts[1]

# test add_no_redundant
# test_mol_fails_roundtrip
# Test file does not exist