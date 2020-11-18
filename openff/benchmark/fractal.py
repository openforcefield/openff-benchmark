"""
Top-level QCFractal cluster management components.

These are largely wrappers around the control points of QCFractal,
with many configuration details either hidden or hard-set for our use-case.

"""

from qcfractal.cli.qcfractal_server import server_init, server_start


def fractal_server_init():
    """Initialize QCFractal server, including any configuration items specific to benchmarking.

    """
    # be sure to initialize DB in a reliable place on filesystem
    # default is `~/.qca`, which should be good


def fractal_server_start():
    """Start QCFractal server, including any configuration items specific to benchmarking.

    """
    # be sure to disable authentication


def fractal_manager_start():
    """Start QCFractal manager, including any configuration items specific to benchmarking.

    """
    # need configuration generation
    # will require survey results to get this right
