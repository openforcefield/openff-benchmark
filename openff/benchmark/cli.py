"""
Top-level `openff-benchmark` cli entrypoint.

"""

import click

                    


@click.group()
def cli():
    pass

@cli.group()
def fractal():
    """Control points for a fractal server and managers.

    """
    pass

@fractal.command()
def server_init():
    from .fractal import fractal_server_init
    fractal_server_init()

@fractal.command()
def server_start():
    from .fractal import fractal_server_start
    fractal_server_start()

@fractal.command()
def manager_init():
    pass

@fractal.command()
def manager_start():
    from .fractal import fractal_manager_start
    fractal_manager_start()

@cli.group()
def optimize():
    """Execute benchmarking geometry optimizations.

    """
    pass

@optimize.command()
@click.option('--server-uri', default="localhost:7777")
@click.option('--dataset-name', default="Benchmark Optimizations")
@click.option('--season', default=1)
@click.argument('input-path', nargs=-1)
def submit_molecules(server_uri, input_path, season, dataset_name):
    from .geometry_optimizations import compute as optcompute
    optcompute.submit_molecules(
            server_uri, input_path, season, dataset_name)


@optimize.command()
def submit_compute(server_uri, input_path):
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
