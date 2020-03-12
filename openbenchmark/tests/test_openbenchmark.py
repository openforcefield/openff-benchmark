"""
Unit and regression test for the openbenchmark package.
"""

# Import package, test suite, and other packages as needed
import openbenchmark
import pytest
import sys

def test_openbenchmark_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "openbenchmark" in sys.modules
