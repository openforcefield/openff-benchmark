from openforcefield.utils.toolkits import GLOBAL_TOOLKIT_REGISTRY, OpenEyeToolkitWrapper, RDKitToolkitWrapper
from openff.benchmark.utils.validate_and_assign_ids import validate_and_assign
from openff.benchmark.utils.utils import get_data_file_path
import pytest

def test_openeye_deregistered():
    for toolkitwrapper in GLOBAL_TOOLKIT_REGISTRY.registered_toolkits:
        assert not isinstance(toolkitwrapper, OpenEyeToolkitWrapper)
    to_smiles = GLOBAL_TOOLKIT_REGISTRY.resolve('to_smiles')
    assert not isinstance(to_smiles.__self__, OpenEyeToolkitWrapper)
    assert isinstance(to_smiles.__self__, RDKitToolkitWrapper)
#    raise Exception(dir(to_smiles))
        
def test_validate_input_good():
    validate_and_assign(get_data_file_path('input_good.sdf'), 'aaa')
            
def test_validate_input_single_mol():
    validate_and_assign(get_data_file_path('input_single_mol.sdf'), 'aaa')
        
def test_validate_input_stereoisomers():
    validate_and_assign(get_data_file_path('input_stereoisomers.sdf'), 'aaa')
    
def test_validate_input_duplicates():
    with pytest.raises(Exception, match="Duplicate"):
        validate_and_assign(get_data_file_path('input_duplicates.sdf'), 'aaa')

def test_validate_input_2d():
    with pytest.raises(Exception, match="stereochemistry"):
        validate_and_assign(get_data_file_path('input_2d.sdf'), 'aaa')



        
