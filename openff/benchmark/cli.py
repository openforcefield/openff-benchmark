"""
Top-level `openff-benchmark` cli entrypoint.

"""

import click

from .fractal import (fractal_server_init,
                      fractal_server_start,
                      fractal_manager_start)

from .geometry_optimizations import compute as optcompute
                    


@click.group()
def cli():
    pass

@cli.group()
def fractal():
    pass

@fractal.command()
def server_init():
    fractal_server_init()

@fractal.command()
def server_start():
    fractal_server_start()

@fractal.command()
def manager_init():
    pass

@fractal.command()
def manager_start():
    pass

@cli.group()
def optimize():
    """Control points for executing benchmarking geometry optimizations.

    """
    pass

@optimize.command()
def submit_molecules():
    pass

@optimize.command()
def submit_compute():
    pass

@optimize.command()
def export():
    pass

@optimize.command()
def status():
    pass

@optimize.command()
@click.option('--mode',
              type=click.Choice(['singleshot', 'service']),
              default='singleshot',
              help="If 'singleshot', runs once and exits; if 'service', runs in a loop")
@click.option('--sleep-time', default=3600, help="Time in seconds to sleep when in 'service' mode")
@click.option('--only-compute', default=None, help="Comma-separated compute spec names to limit errorcycling to")
@click.argument('molids', nargs=-1)
def errorcycle(molids, mode, sleep_time, only_compute):
    pass

@optimize.command()
def local_from_server():
    pass

@optimize.command()
def local_from_molecule():
    pass
