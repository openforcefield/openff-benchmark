from openforcefield.utils.toolkits import GLOBAL_TOOLKIT_REGISTRY, OpenEyeToolkitWrapper, RDKitToolkitWrapper
from openff.benchmark.utils.validate_and_assign_ids import validate_and_assign
from openff.benchmark.utils.utils import get_data_file_path, temporary_cd
import pytest
import tempfile
import os
import shutil
import inspect
import glob

def test_openeye_deregistered():
    for toolkitwrapper in GLOBAL_TOOLKIT_REGISTRY.registered_toolkits:
        assert not isinstance(toolkitwrapper, OpenEyeToolkitWrapper)
    to_smiles = GLOBAL_TOOLKIT_REGISTRY.resolve('to_smiles')
    assert not isinstance(to_smiles.__self__, OpenEyeToolkitWrapper)
    assert isinstance(to_smiles.__self__, RDKitToolkitWrapper)

def test_dont_overwrite_output_directory(tmpdir):
    with tmpdir.as_cwd():
        test_name = inspect.stack()[0].function
        test_dir = os.path.join(test_name, '1-validate_and_assign_graphs_and_confs')
        if os.path.exists(test_dir):
            shutil.rmtree(test_dir)
        input_mols = [get_data_file_path('input_single_mol.sdf')]
        validate_and_assign(input_mols,
                            'BBB',
                            test_dir,
                            )
        with pytest.raises(Exception, match='Specify `delete_existing=True` to remove'):
            validate_and_assign(input_mols,
                                'BBB',
                                test_dir
                                )


# Conformers w/o graphs
class Test3dInputsWOGraphs:
    # single file single mol
    def test_single_file_single_mol(self, tmpdir):
        with tmpdir.as_cwd():
            test_name = f'{self.__class__.__name__}.{inspect.stack()[0].function}'
            test_dir = os.path.join(test_name, '1-validate_and_assign_graphs_and_confs')
            #test_dir = 'test_multi_molecule_multi_sdf_graph_input'
            if os.path.exists(test_dir):
                shutil.rmtree(test_dir)
            #os.makedirs(test_dir)
            input_mols = [get_data_file_path('input_single_mol.sdf')]
            validate_and_assign(input_mols,
                                'BBB',
                                test_dir,
                                )
            output_files = glob.glob(os.path.join(test_dir, '*.sdf'))
            output_files = [os.path.basename(fname) for fname in output_files]
            #assert 'BBB-00000.smi' in output_files
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
            test_name = f'{self.__class__.__name__}.{inspect.stack()[0].function}'
            test_dir = os.path.join(test_name, '1-validate_and_assign_graphs_and_confs')
            #test_dir = 'test_multi_molecule_multi_sdf_graph_input'
            if os.path.exists(test_dir):
                shutil.rmtree(test_dir)
            #os.makedirs(test_dir)
            input_mols = [get_data_file_path('input_multi_mol.sdf')]
            validate_and_assign(input_mols,
                                'BBB',
                                test_dir,
                                )
            output_files = glob.glob(os.path.join(test_dir, '*.sdf'))
            output_files = [os.path.basename(fname) for fname in output_files]
            #assert 'BBB-00000.smi' in output_files
            #assert 'BBB-00003.smi' in output_files
            assert 'BBB-00000-00.sdf' in output_files
            assert 'BBB-00003-00.sdf' in output_files
            assert len(output_files) == 4
        
    # single file multi conformer
    def test_single_file_multi_conformer(self, tmpdir):
        with tmpdir.as_cwd():
            test_name = f'{self.__class__.__name__}.{inspect.stack()[0].function}'
            test_dir = os.path.join(test_name, '1-validate_and_assign_graphs_and_confs')
            #test_dir = 'test_multi_molecule_multi_sdf_graph_input'
            if os.path.exists(test_dir):
                shutil.rmtree(test_dir)
            #os.makedirs(test_dir)
            input_mols = [get_data_file_path('input_duplicates.sdf')]
            validate_and_assign(input_mols,
                                'BBB',
                                test_dir,
                                )
            output_files = glob.glob(os.path.join(test_dir, '*.sdf'))
            output_files = [os.path.basename(fname) for fname in output_files]
            #assert 'BBB-00000.smi' in output_files
            assert 'BBB-00000-00.sdf' in output_files
            assert 'BBB-00000-09.sdf' in output_files
            assert len(output_files) == 10
        
    # multi file single mol
    @pytest.mark.xfail
    def test_multi_file_single_mol(self, tmpdir):
        with tmpdir.as_cwd():
            test_name = f'{self.__class__.__name__}.{inspect.stack()[0].function}'
            test_dir = os.path.join(test_name, '1-validate_and_assign_graphs_and_confs')
            #test_dir = 'test_multi_molecule_multi_sdf_graph_input'
            if os.path.exists(test_dir):
                shutil.rmtree(test_dir)
            #os.makedirs(test_dir)
            input_mols = [get_data_file_path('input_single_mol.sdf'),
                          get_data_file_path('input_one_stereoisomer.sdf')]
            validate_and_assign(input_mols,
                                'BBB',
                                test_dir,
                                )
            output_files = glob.glob(os.path.join(test_dir, '*.sdf'))
            output_files = [os.path.basename(fname) for fname in output_files]
            #assert 'BBB-00000.smi' in output_files
            #assert 'BBB-00001.smi' in output_files
            assert 'BBB-00000-00.sdf' in output_files
            assert 'BBB-00001-00.sdf' in output_files
            assert len(output_files) == 2
        
    # multi file multi mol
    @pytest.mark.xfail
    def test_multi_file_multi_mol(self, tmpdir):
        with tmpdir.as_cwd():
            test_name = f'{self.__class__.__name__}.{inspect.stack()[0].function}'
            test_dir = os.path.join(test_name, '1-validate_and_assign_graphs_and_confs')
            #test_dir = 'test_multi_molecule_multi_sdf_graph_input'
            if os.path.exists(test_dir):
                shutil.rmtree(test_dir)
            #os.makedirs(test_dir)
            input_mols = [get_data_file_path('input_multi_mol.sdf'),
                          get_data_file_path('input_all_stereoisomers.sdf')]
            validate_and_assign(input_mols,
                                'BBB',
                                test_dir,
                                )
            output_files = glob.glob(os.path.join(test_dir, '*.sdf'))
            output_files = [os.path.basename(fname) for fname in output_files]
            #assert 'BBB-00000.smi' in output_files
            #assert 'BBB-00011.smi' in output_files
            #assert 'BBB-00012.smi' not in output_files
            assert 'BBB-00000-00.sdf' in output_files
            assert 'BBB-00011-00.sdf' in output_files
            assert len(output_files) == 12

    # multi file multi mol duplicates (ok)
    pytest.mark.xfail    
    def test_multi_file_multi_mol_duplicates(self, tmpdir):
        with tmpdir.as_cwd():
            #raise Exception(inspect.stack()[0])
            #test_name = f'{self.__class__.__name__}.{inspect.stack()[0].function}'
            test_name = f'{self.__class__.__name__}.{inspect.stack()[0].function}'
            test_dir = os.path.join(test_name, '1-validate_and_assign_graphs_and_confs')
            #test_dir = 'test_multi_molecule_multi_sdf_graph_input'
            if os.path.exists(test_dir):
                shutil.rmtree(test_dir)
            #os.makedirs(test_dir)
            input_mols = [get_data_file_path('input_one_stereoisomer.sdf'),
                          get_data_file_path('input_all_stereoisomers.sdf')]
            validate_and_assign(input_mols,
                                'BBB',
                                test_dir,
                                )
            output_files = glob.glob(os.path.join(test_dir, '*.sdf'))
            output_files = [os.path.basename(fname) for fname in output_files]
            #assert 'BBB-00000.smi' in output_files
            #assert 'BBB-00007.smi' in output_files
            assert 'BBB-00000-00.sdf' in output_files
            assert 'BBB-00000-01.sdf' not in output_files
            assert 'BBB-00007-00.sdf' in output_files
            assert 'BBB-00008-00.sdf' not in output_files
            assert len(output_files) == 8

#
# Conformers w graphs
# class TestGraphAnd3dInputs:
#     # no overlap
#     def test_no_overlap(self, tmpdir):
#         with tmpdir.as_cwd():
#             test_name = f'{self.__class__.__name__}.{inspect.stack()[0].function}'
#             test_dir = os.path.join(test_name, '1-validate_and_assign_graphs_and_confs')
#             #test_dir = 'test_multi_molecule_multi_sdf_graph_input'
#             if os.path.exists(test_dir):
#                 shutil.rmtree(test_dir)
#             #os.makedirs(test_dir)
#             input_graphs = [get_data_file_path('input_single_mol.sdf')]
#             input_3ds = [get_data_file_path('input_one_stereoisomer.sdf')]
#             validate_and_assign(input_graphs,
#                                 input_3ds,
#                                 'BBB',
#                                 test_dir,
#                                 )
#             output_files = glob.glob(os.path.join(test_dir, '*'))
#             output_files = [os.path.basename(fname) for fname in output_files]
#             assert 'BBB-00000.smi' in output_files
#             assert 'BBB-00000-00.sdf' not in output_files
#             assert 'BBB-00001.smi' in output_files
#             assert 'BBB-00001-00.sdf' in output_files
#             assert len(output_files) == 3
        
#     # graphs superset of conformers
#     def test_graphs_superset_of_conformers(self, tmpdir):
#         with tmpdir.as_cwd():
#             test_name = f'{self.__class__.__name__}.{inspect.stack()[0].function}'
#             test_dir = os.path.join(test_name, '1-validate_and_assign_graphs_and_confs')
#             #test_dir = 'test_multi_molecule_multi_sdf_graph_input'
#             if os.path.exists(test_dir):
#                 shutil.rmtree(test_dir)
#             #os.makedirs(test_dir)
#             input_graphs = [get_data_file_path('input_single_mol.sdf'),
#                             get_data_file_path('input_one_stereoisomer.sdf')]
#             input_3ds = [get_data_file_path('input_one_stereoisomer.sdf')]
#             validate_and_assign(input_graphs,
#                                 input_3ds,
#                                 'BBB',
#                                 test_dir,
#                                 )
#             output_files = glob.glob(os.path.join(test_dir, '*'))
#             output_files = [os.path.basename(fname) for fname in output_files]
#             assert 'BBB-00000.smi' in output_files
#             assert 'BBB-00000-00.sdf' not in output_files
#             assert 'BBB-00001.smi' in output_files
#             assert 'BBB-00001-00.sdf' in output_files
#             assert len(output_files) == 3
        
#     # conformers superset of graphs
#     def test_conformers_superset_of_graphs(self, tmpdir):
#         with tmpdir.as_cwd():
#             test_name = f'{self.__class__.__name__}.{inspect.stack()[0].function}'
#             test_dir = os.path.join(test_name, '1-validate_and_assign_graphs_and_confs')
#             #test_dir = 'test_multi_molecule_multi_sdf_graph_input'
#             if os.path.exists(test_dir):
#                 shutil.rmtree(test_dir)
#             #os.makedirs(test_dir)
#             input_graphs = [get_data_file_path('input_single_mol.sdf')]
#             input_3ds = [get_data_file_path('input_one_stereoisomer.sdf'),
#                          get_data_file_path('input_single_mol.sdf')]
#             validate_and_assign(input_graphs,
#                                 input_3ds,
#                                 'BBB',
#                                 test_dir,
#                                 )
#             output_files = glob.glob(os.path.join(test_dir, '*'))
#             output_files = [os.path.basename(fname) for fname in output_files]
#             assert 'BBB-00000.smi' in output_files
#             assert 'BBB-00000-00.sdf' in output_files
#             assert 'BBB-00001.smi' in output_files
#             assert 'BBB-00001-00.sdf' in output_files
#             assert len(output_files) == 4
        
