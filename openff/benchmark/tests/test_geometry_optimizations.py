"""Tests for geometry optimization library components.

"""

import pytest

from openff.benchmark.geometry_optimizations import compute
from openff.benchmark.utils.utils import get_data_file_path

def test_submit_molecules():
    pass

def test_export_molecule_data():
    pass

def test_get_optimization_status():
    pass

def test_errorcycle_optimizations():
    pass

def test_optimization_from_server():
    pass

def test_execute_optimization_from_molecules(tmpdir):
    with tmpdir.as_cwd():

        input_mols = [get_data_file_path('ethane.sdf')]

        compute.execute_optimization_from_molecules(
                input_mols, 
                '3-export-compute',
                '1:1')
