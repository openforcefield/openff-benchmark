from openff.benchmark.utils.generate_conformers import generate_conformers
import inspect
import os
from openff.benchmark.utils.utils import get_data_file_path

# test loading a mix of graph and 3D molecules, generating up to 10 confs
def test_generate_conformers():
    test_name = inspect.stack()[0].function
    test_dir = os.path.join(test_name, '2-generate_conformers') 
    generate_conformers(get_data_file_path('1-validate_and_assign'),
                        test_dir)

# test loading a molecule with 10 pre-set conformers, so no new confs should be generated

# test loading a molecule with pre-set conformers within the RMSD threshold of each other


