"""
Unit and regression test for the openbenchmark package.
"""

# Import package, test suite, and other packages as needed
from openff import benchmark
import pytest
import sys

def test_openff_benchmark_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "openff.benchmark" in sys.modules
