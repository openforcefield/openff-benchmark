from openff.benchmark.utils.coverage_report import generate_coverage_report
from openff.benchmark.utils.utils import get_data_file_path
import pytest

# Test single file
@pytest.mark.parametrize('n_procs', [1,2])
def test_single_file(n_procs):
    in_file = get_data_file_path('1-validate_and_assign_graphs_and_confs/BBB-00000-00.sdf')
    coverage, success_mols, error_mols = generate_coverage_report(input_molecules=in_file,
                                                                  forcefield_name='openff_unconstrained-1.3.0.offxml',
                                                                  processors=n_procs)
    assert len(success_mols) == 1, error_mols[0][1]
    assert len(error_mols) == 0

# Test multiple files
@pytest.mark.parametrize('n_procs', [1,2])
def test_multiple_files(n_procs):
    in_files = [get_data_file_path('1-validate_and_assign_graphs_and_confs/BBB-00000-00.sdf'),
               get_data_file_path('1-validate_and_assign_graphs_and_confs/BBB-00001-00.sdf')]
    coverage, success_mols, error_mols = generate_coverage_report(input_molecules=in_files,
                                                                  forcefield_name='openff_unconstrained-1.3.0.offxml',
                                                                  processors=n_procs)
    assert len(success_mols) == 2
    assert len(error_mols) == 0
    #raise Exception(error_mols)
# Test not double counting if multiple conformers in the same directory
def test_input_directory():
    in_files = get_data_file_path('1-validate_and_assign_graphs_and_confs/')
    coverage, success_mols, error_mols = generate_coverage_report(input_molecules=in_files,
                                                                  forcefield_name='openff_unconstrained-1.3.0.offxml')
    assert len(success_mols) == 4
    assert len(error_mols) == 0

# Test catching uncovered params
@pytest.mark.parametrize('n_procs',[1,2])
def test_error_missing_valence_param(n_procs):
    in_file = get_data_file_path('missing_valence_params.sdf')
    coverage, success_mols, error_mols = generate_coverage_report(input_molecules=in_file,
                                                                  forcefield_name='openff_unconstrained-1.3.0.offxml',
                                                                  processors=n_procs)
    assert len(success_mols) == 0
    assert len(error_mols) == 1
    assert "BondHandler was not able to find parameters" in str(error_mols[0][1])

# Test for a molecule covered by OpenFF but that causes antechamber to fail
# I can't actually come up with such a molecule -- Suggestions welcome
# Also note that charge gen may be getting skipped --
# Search "deregister_parameter_handler" in coverage_report.py
@pytest.mark.skip
def test_error_uncovered_antechamber_param():
    in_file = get_data_file_path('sodium_carbide.sdf')
    coverage, success_mols, error_mols = generate_coverage_report(input_molecules=in_file,
                                                                  forcefield_name='openff_unconstrained-1.3.0.offxml')
    assert len(success_mols) == 0
    assert len(error_mols) == 1
    assert "Command '['antechamber'" in str(error_mols[0][1])
