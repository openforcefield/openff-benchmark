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
#    raise Exception(dir(to_smiles))
        
# def test_validate_input_good():
#     validate_and_assign(get_data_file_path('input_good.sdf'), 'aaa')
            
# def test_validate_input_single_mol():
#     validate_and_assign(get_data_file_path('input_single_mol.sdf'), 'aaa')
        
# def test_validate_input_stereoisomers():
#     validate_and_assign(get_data_file_path('input_stereoisomers.sdf'), 'aaa')
    
# def test_validate_input_duplicates():
#     with pytest.raises(Exception, match="Duplicate"):
#         validate_and_assign(get_data_file_path('input_duplicates.sdf'), 'aaa')

# def test_validate_input_2d():
#     with pytest.raises(Exception, match="stereochemistry"):
#         validate_and_assign(get_data_file_path('input_2d.sdf'), 'aaa')

# def test_single_molecule_single_sdf_graph_input():
#     test_dir = 'test_single_molecule_graph_input'
#     if os.path.exists(test_dir):
#         shutil.rmtree(test_dir)
#     os.mkdir(test_dir)

#     input_mols = [get_data_file_path('input_single.sdf')]
#     validate_and_assign(input_mols,
#                         '',
#                         os.path.join(test_dir, '1-validate_and_assign'),
#                         'BBB')

# def test_multi_molecule_single_sdf_graph_input():
#     test_dir = 'test_single_molecule_graph_input'
#     if os.path.exists(test_dir):
#         shutil.rmtree(test_dir)
#     os.mkdir(test_dir)

#     input_mols = [get_data_file_path('input_multimolecule.sdf')]
#     validate_and_assign(input_mols,
#                         '',
#                         os.path.join(test_dir, '1-validate_and_assign'),
#                         'BBB')
    
# def test_multi_file_multi_molecule_graph_input():
#     #with tempfile.TemporaryDirectory() as tmpdir:
#     #    with temporary_cd(tmpdir):
#     test_dir = 'test_multi_molecule_multi_sdf_graph_input'
#     if os.path.exists(test_dir):
#         shutil.rmtree(test_dir)
#     os.mkdir(test_dir)
#     input_mols = [get_data_file_path('input_good.sdf'),
#                   get_data_file_path('input_stereoisomers.sdf')]
#     validate_and_assign(input_mols,
#                         '',
#                         os.path.join(test_dir, '1-validate_and_assign'),
#                         'AAA')
    
# def test_multi_molecule_i():
#     #with tempfile.TemporaryDirectory() as tmpdir:
#     #    with temporary_cd(tmpdir):
#     test_dir = 'test_multi_molecule_multi_sdf_graph_input'
#     if os.path.exists(test_dir):
#         shutil.rmtree(test_dir)
#     os.mkdir(test_dir)
#     input_mols = [get_data_file_path('input_good.sdf'),
#                   get_data_file_path('input_stereoisomers.sdf')]
#     validate_and_assign(input_mols,
#                         '',
#                         os.path.join(test_dir, '1-validate_and_assign'),
#                         'AAA')


#@contextlib.contextmanager
#def clean_test_subdir():
#    test_name = inspect.stack()[1][3]
#    #prev_dir = os.getcwd()
#    #os.chdir(os.path.abspath(dir_path))
#    if os.path.exists(test_dir):
#        shutil.rmtree(test_dir)
#    os.mkdir(test_dir)
#    yield
#    #try:
#    #    yield
#    #finally:
#    #    os.chdir(prev_dir)

# Graph inputs w/o conformers
class TestGraphInputsWOConformers:
    # single file single mol
    def test_single_file_single_mol(self):
        #test_name = inspect.stack()[1][3]
        test_name = inspect.stack()[0].function
        test_dir = os.path.join(test_name, '1-validate_and_assign')
        #test_dir = 'test_multi_molecule_multi_sdf_graph_input'
        if os.path.exists(test_dir):
            shutil.rmtree(test_dir)
        #os.makedirs(test_dir)
        input_mols = [get_data_file_path('input_single_mol.sdf')]
        validate_and_assign(input_mols,
                            '',
                            test_dir,
                            'BBB')
        output_files = glob.glob(os.path.join(test_dir, '*'))
        output_files = [os.path.basename(fname) for fname in output_files]
        assert 'BBB-00000.smi' in output_files
        assert len(output_files) == 1
        
    # single file multi mol
    def test_single_file_multi_mol(self):
        test_name = inspect.stack()[0].function
        #raise Exception(inspect.stack()[0].function)
        test_dir = os.path.join(test_name, '1-validate_and_assign')
        #test_dir = 'test_multi_molecule_multi_sdf_graph_input'
        if os.path.exists(test_dir):
            shutil.rmtree(test_dir)
        #os.makedirs(test_dir)
        input_mols = [get_data_file_path('input_multi_mol.sdf')]
        validate_and_assign(input_mols,
                            '',
                            test_dir,
                            'BBB')
        output_files = glob.glob(os.path.join(test_dir, '*'))
        output_files = [os.path.basename(fname) for fname in output_files]
        assert 'BBB-00000.smi' in output_files
        assert 'BBB-00003.smi' in output_files
        assert len(output_files) == 4
        
    # multi file single mol
    def test_multi_file_single_mol(self):
        test_name = inspect.stack()[0].function
        #raise Exception(inspect.stack()[0].function)
        test_dir = os.path.join(test_name, '1-validate_and_assign')
        #test_dir = 'test_multi_molecule_multi_sdf_graph_input'
        if os.path.exists(test_dir):
            shutil.rmtree(test_dir)
        #os.makedirs(test_dir)
        input_mols = [get_data_file_path('input_single_mol.sdf'),
                      get_data_file_path('input_one_stereoisomer.sdf')]
        validate_and_assign(input_mols,
                            '',
                            test_dir,
                            'BBB')
        output_files = glob.glob(os.path.join(test_dir, '*'))
        output_files = [os.path.basename(fname) for fname in output_files]
        assert 'BBB-00000.smi' in output_files
        assert 'BBB-00001.smi' in output_files
        assert len(output_files) == 2
        
    # multi file multi mol
    def test_multi_file_multi_mol(self):
        test_name = inspect.stack()[0].function
        #raise Exception(inspect.stack()[0].function)
        test_dir = os.path.join(test_name, '1-validate_and_assign')
        #test_dir = 'test_multi_molecule_multi_sdf_graph_input'
        if os.path.exists(test_dir):
            shutil.rmtree(test_dir)
        #os.makedirs(test_dir)
        input_mols = [get_data_file_path('input_multi_mol.sdf'),
                      get_data_file_path('input_all_stereoisomers.sdf')]
        validate_and_assign(input_mols,
                            '',
                            test_dir,
                            'BBB')
        output_files = glob.glob(os.path.join(test_dir, '*'))
        output_files = [os.path.basename(fname) for fname in output_files]
        assert 'BBB-00000.smi' in output_files
        assert 'BBB-00011.smi' in output_files
        assert len(output_files) == 12
        
    # single file multi mol duplicates (error)
    def test_single_file_multi_mol_duplicate_error(self):
        test_name = inspect.stack()[0].function
        #raise Exception(inspect.stack()[0].function)
        test_dir = os.path.join(test_name, '1-validate_and_assign')
        #test_dir = 'test_multi_molecule_multi_sdf_graph_input'
        if os.path.exists(test_dir):
            shutil.rmtree(test_dir)
        #os.makedirs(test_dir)
        input_mols = [get_data_file_path('input_duplicates.sdf')]
        with pytest.raises(Exception, match='Duplicate'):
            validate_and_assign(input_mols,
                                '',
                                test_dir,
                                'BBB')
            
    # multi file multi mol duplicates (error)
    def test_multi_file_multi_mol_duplicate_error(self):
        test_name = inspect.stack()[0].function
        #raise Exception(inspect.stack()[0].function)
        test_dir = os.path.join(test_name, '1-validate_and_assign')
        #test_dir = 'test_multi_molecule_multi_sdf_graph_input'
        if os.path.exists(test_dir):
            shutil.rmtree(test_dir)
        #os.makedirs(test_dir)
        input_mols = [get_data_file_path('input_one_stereoisomer.sdf'),
                      get_data_file_path('input_all_stereoisomers.sdf')]
        with pytest.raises(Exception, match='Duplicate'):
            validate_and_assign(input_mols,
                                '',
                                test_dir,
                                'BBB')
            
    # input undefined stereochemistry (error)
    def test_undefined_stereochemistry(self):
        test_name = inspect.stack()[0].function
        #raise Exception(inspect.stack()[0].function)
        test_dir = os.path.join(test_name, '1-validate_and_assign')
        #test_dir = 'test_multi_molecule_multi_sdf_graph_input'
        if os.path.exists(test_dir):
            shutil.rmtree(test_dir)
        #os.makedirs(test_dir)
        input_mols = [get_data_file_path('input_2d.sdf')]
        with pytest.raises(Exception, match='stereo'):
            validate_and_assign(input_mols,
                                '',
                                test_dir,
                                'BBB')

# Conformers w/o graphs
class Test3dInputsWOGraphs:
    # single file single mol
    def test_single_file_single_mol(self):
        test_name = inspect.stack()[0].function
        test_dir = os.path.join(test_name, '1-validate_and_assign')
        #test_dir = 'test_multi_molecule_multi_sdf_graph_input'
        if os.path.exists(test_dir):
            shutil.rmtree(test_dir)
        #os.makedirs(test_dir)
        input_mols = [get_data_file_path('input_single_mol.sdf')]
        validate_and_assign('',
                            input_mols,
                            test_dir,
                            'BBB')
        output_files = glob.glob(os.path.join(test_dir, '*'))
        output_files = [os.path.basename(fname) for fname in output_files]
        assert 'BBB-00000.smi' in output_files
        assert 'BBB-00000-000.sdf' in output_files
        assert len(output_files) == 2
    
    # single file multi mol
    def test_single_file_multi_mol(self):
        test_name = inspect.stack()[0].function
        test_dir = os.path.join(test_name, '1-validate_and_assign')
        #test_dir = 'test_multi_molecule_multi_sdf_graph_input'
        if os.path.exists(test_dir):
            shutil.rmtree(test_dir)
        #os.makedirs(test_dir)
        input_mols = [get_data_file_path('input_multi_mol.sdf')]
        validate_and_assign('',
                            input_mols,
                            test_dir,
                            'BBB')
        output_files = glob.glob(os.path.join(test_dir, '*'))
        output_files = [os.path.basename(fname) for fname in output_files]
        assert 'BBB-00000.smi' in output_files
        assert 'BBB-00000-000.sdf' in output_files
        assert 'BBB-00000-003.sdf' in output_files
        assert len(output_files) == 5
        
#     multi file single mol
#     multi file multi mol
#     single file multi mol duplicates (ok)
#     multi file multi mol duplicates (ok)
#
# Conformers w graphs
#    graphs superset of conformers
#    conformers superset of graphs
    
