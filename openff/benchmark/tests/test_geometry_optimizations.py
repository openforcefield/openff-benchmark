"""Tests for geometry optimization library components.

"""
import os
import re

import pytest
import pint
from click.testing import CliRunner

from openforcefield.topology import Molecule

from openff.benchmark.geometry_optimizations import compute
from openff.benchmark.geometry_optimizations.seasons import SEASONS
from openff.benchmark.utils.utils import get_data_file_path
from openff.benchmark.cli import cli

from qcportal import FractalClient

runner = CliRunner()

punit = pint.UnitRegistry()
Q_ = punit.Quantity

# CLI tests

@pytest.fixture(scope="module")
def ethane_preprocessed():
    """Ethane, with prepopulated properties for identification.

    The `group_name`, `molecule_index`, and conformer_index` fields are set as:

    ```
    >  <group_name>  (1) 
    OFF
    
    >  <molecule_index>  (1) 
    00000
    
    >  <conformer_index>  (1) 
    00
    ```

    """
    return get_data_file_path('ethane.sdf')


@pytest.fixture(scope="module")
def ethane_qm_optimized():
    """Ethane, QM optimized.

    The `group_name`, `molecule_index`, and conformer_index` fields are set as:

    ```
    >  <group_name>  (1) 
    OFF
    
    >  <molecule_index>  (1) 
    00000
    
    >  <conformer_index>  (1) 
    00
    ```

    """
    return get_data_file_path('ethane-qm-optimized.sdf')


@pytest.fixture(scope="module")
def acetic_acid_optimized_proton_transfer():
    """Acetic acid, QM optimized, with the connection table indicating a proton transfer between oxygens.

    The `group_name`, `molecule_index`, and conformer_index` fields are set as:

    ```
    >  <group_name>  (1)
    BBB

    >  <molecule_index>  (1)
    00000

    >  <conformer_index>  (1)
    00
    ```

    """
    return get_data_file_path('acetic_acid_optimized_proton_transfer.sdf')

## server required

def test_cli_optimize_submit_molecules_qm(fractal_compute_server, tmpdir, ethane_preprocessed):
    """Test submission of QM to a Fractal server via the CLI.

    """

    fc = FractalClient(fractal_compute_server)

    with tmpdir.as_cwd():
        response = runner.invoke(cli, ['optimize', 'submit-molecules',
                                       '--fractal-uri', fc.address,
                                       '--dataset-name', '"Test Dataset - QM"',
                                       '--season', '1:1',
                                       ethane_preprocessed])

        assert response.exit_code == 0
        assert "Submitted!" in response.output
    
        # wait for results
        fractal_compute_server.await_results()

        response = runner.invoke(cli, ['optimize', 'status',
                                       '--fractal-uri', fc.address,
                                       '--dataset-name', '"Test Dataset - QM"'])

        assert response.exit_code == 0
        assert "COMPLETE" in response.output


def test_cli_optimize_submit_molecules_mm(fractal_compute_server, tmpdir, ethane_qm_optimized):
    """Test submission of QM to a Fractal server via the CLI.

    """

    fc = FractalClient(fractal_compute_server)

    with tmpdir.as_cwd():
        response = runner.invoke(cli, ['optimize', 'submit-molecules',
                                       '--fractal-uri', fc.address,
                                       '--dataset-name', '"Test Dataset - MM"',
                                       '--season', '1:2',
                                       ethane_qm_optimized])

        assert response.exit_code == 0
        assert "Submitted!" in response.output
    
        # wait for results
        fractal_compute_server.await_results()

        response = runner.invoke(cli, ['optimize', 'status',
                                       '--fractal-uri', fc.address,
                                       '--dataset-name', '"Test Dataset - MM"'])

        assert response.exit_code == 0
        assert "COMPLETE" in response.output


def test_cli_optimize_export_qm(fractal_compute_server, tmpdir, ethane_preprocessed):
    """Test export of molecule QM data from a Fractal server via the CLI.

    """

    # need to submit first; should be quick if already there
    test_cli_optimize_submit_molecules_qm(fractal_compute_server, tmpdir, ethane_preprocessed)

    fc = FractalClient(fractal_compute_server)

    with tmpdir.as_cwd():
        outdir = '4-compute-qm'
        response = runner.invoke(cli, ['optimize', 'export',
                                       '--fractal-uri', fc.address,
                                       '--dataset-name', '"Test Dataset - QM"',
                                       '-o', outdir])

        assert response.exit_code == 0
        assert "exporting COMPLETE" in response.output

        for season in SEASONS["1:1"]:
            exported_files = os.listdir(os.path.join(tmpdir.strpath, outdir, season))
            assert set(exported_files) - set(['error_mols']) == set(f"OFF-00000-00.{suffix}" for suffix in ["sdf", "json", "perf.json"])


def test_cli_optimize_export_mm(fractal_compute_server, tmpdir, ethane_qm_optimized):
    """Test export of molecule MM data from a Fractal server via the CLI.

    """

    # need to submit first; should be quick if already there
    test_cli_optimize_submit_molecules_mm(fractal_compute_server, tmpdir, ethane_qm_optimized)

    fc = FractalClient(fractal_compute_server)

    with tmpdir.as_cwd():
        outdir = '4-compute-qm'
        response = runner.invoke(cli, ['optimize', 'export',
                                       '--fractal-uri', fc.address,
                                       '--dataset-name', '"Test Dataset - MM"',
                                       '-o', outdir])

        assert response.exit_code == 0
        assert "exporting COMPLETE" in response.output

        for season in SEASONS["1:2"]:
            exported_files = os.listdir(os.path.join(tmpdir.strpath, outdir, season))
            assert set(exported_files) - set(['error_mols']) == set(f"OFF-00000-00.{suffix}" for suffix in ["sdf", "json", "perf.json"])


def test_cli_optimize_status(fractal_compute_server, tmpdir, ethane_preprocessed):
    """Test status output from a Fractal server via the CLI.

    """

    # need to submit first; should be quick if already there
    test_cli_optimize_submit_molecules_qm(fractal_compute_server, tmpdir, ethane_preprocessed)

    fc = FractalClient(fractal_compute_server)

    with tmpdir.as_cwd():
        response = runner.invoke(cli, ['optimize', 'status',
                                       '--fractal-uri', fc.address,
                                       '--dataset-name', '"Test Dataset - QM"'])

        assert response.exit_code == 0
        assert re.compile("OFF-00000-00(\s*)COMPLETE").match(response.output.strip().split('\n')[-1])


def test_cli_optimize_progress():
    pass


## no server required

def test_cli_optimize_execute_qm(tmpdir, ethane_preprocessed):
    """Test execution of QM via QCEngine via the CLI.

    """

    with tmpdir.as_cwd():
        outdir = '4-compute-qm'
        response = runner.invoke(cli, ['optimize', 'execute',
                                       '--season', '1:1',
                                       '--nthreads', '2',
                                       '--memory', '1',
                                       '-o', outdir,
                                       ethane_preprocessed])

        assert response.exit_code == 0
        assert "... 'OFF-00000-00'" in response.output

        for season in SEASONS["1:1"]:
            exported_files = os.listdir(os.path.join(tmpdir.strpath, outdir, season))
            assert set(exported_files) == set(f"OFF-00000-00.{suffix}" for suffix in ["sdf", "json", "perf.json"])


def test_cli_optimize_execute_qm_proton_transfer(tmpdir, acetic_acid_optimized_proton_transfer):
    """Test execution of QM via QCEngine via the CLI.

    """

    with tmpdir.as_cwd():
        outdir = '4-compute-qm'
        response = runner.invoke(cli, ['optimize', 'execute',
                                       '--season', '1:1',
                                       '--nthreads', '2',
                                       '--memory', '1',
                                       '-o', outdir,
                                       acetic_acid_optimized_proton_transfer])

        assert response.exit_code == 0
        assert "... 'OFF-00001-00'" in response.output
        assert "... 'OFF-00001-00' : export error" in response.output

        for season in SEASONS["1:1"]:
            exported_files = os.listdir(os.path.join(tmpdir.strpath, outdir, season))
            assert set(exported_files) == set(['error_mols'])
            exported_error_files = os.listdir(os.path.join(tmpdir.strpath, outdir, season, 'error_mols'))
            assert set(exported_error_files) == set(f"OFF-00001-00.{suffix}" for suffix in ["txt", "json", "perf.json"])

def test_cli_optimize_execute_mm(tmpdir, ethane_qm_optimized):
    """Test execution of MM via QCEngine via the CLI.

    """

    with tmpdir.as_cwd():
        outdir = '4-compute-mm'
        response = runner.invoke(cli, ['optimize', 'execute',
                                       '--season', '1:2',
                                       '--nthreads', '2',
                                       '--memory', '1',
                                       '-o', outdir,
                                       ethane_qm_optimized])

        assert response.exit_code == 0
        assert "... 'OFF-00000-00'" in response.output

        for season in SEASONS["1:2"]:
            exported_files = os.listdir(os.path.join(tmpdir.strpath, outdir, season))
            assert set(exported_files) == set(f"OFF-00000-00.{suffix}" for suffix in ["sdf", "json", "perf.json"])


class TestOptimizationExecutor:

    ## server required

    def test_export_molecule_data(self):
        pass

    def test_get_optimization_status(self):
        pass

    def test_errorcycle_optimizations(self):
        pass
    
    def test_optimization_from_server(self):
        pass

    ## no server required
    
    def test_execute_optimization_from_molecules(self, tmpdir):
        with tmpdir.as_cwd():
    
            input_mols = [get_data_file_path('ethane.sdf')]
    
            optexec = compute.OptimizationExecutor()
    
            optexec.execute_optimization_from_molecules(
                    input_mols, 
                    '3-export-compute',
                    '1:1',
                    ncores=1)
    
            file = '3-export-compute/b3lyp-d3bj/dzvp/OFF-00000-00.sdf'
            assert os.path.exists(file)
            
            mol = Molecule.from_file(file) 
    
            assert Q_(mol.properties['initial_energy']) > Q_(mol.properties['final_energy'])
