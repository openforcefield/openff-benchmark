"""
openff-benchmark
Comparison benchmarks between public force fields and Open Force Field Initiative force fields
"""

# Add imports here
from .webapp import *

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
