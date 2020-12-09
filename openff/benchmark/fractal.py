"""
Top-level QCFractal cluster management components.

These are largely wrappers around the control points of QCFractal,
with many configuration details either hidden or hard-set for our use-case.

"""

from qcfractal.cli import qcfractal_server, qcfractal_manager


def fractal_server_init():
    """Initialize QCFractal server, including any configuration items specific to benchmarking.

    """
    # be sure to initialize DB in a reliable place on filesystem
    # default is `~/.qca`, which should be good
    args = {"overwrite_config": True,
            "clear_database": True}

    config = qcfractal_server.FractalConfig()
    
    # might not be terribly safe to use this
    qcfractal_server.server_init(args, config)


def fractal_server_start():
    """Start QCFractal server, including any configuration items specific to benchmarking.

    """
    # be sure to disable authentication; not needed, and just a barrier
    args = {"local_manager": False,
            "disable_ssl": False,
            "tls_key": None,
            "tls_cert": None,
            "start_periodics": True}

    config = qcfractal_server.FractalConfig()

    # might not be terribly safe to use this
    qcfractal_server.server_start(args, config)


def fractal_manager_init():
    """Initialize QCFractal manager, including any configuration items specific to benchmarking.

    """
    # need configuration generation
    # will require survey results to get this right
    # could be done as a wizard-style prompt in the cli

def fractal_manager_start():
    """Start QCFractal manager, including any configuration items specific to benchmarking.

    """
    # TODO: NEEDS WORK
    #settings = qcfractal_manager.ManagerSettings()
    #qcfractal_manager.main()
