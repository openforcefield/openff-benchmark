"""
Top-level `openff-benchmark` cli entrypoint.

"""

import click

from .fractal import (fractal_server_init,
                      fractal_server_start,
                      fractal_manager_start)

from .geometry_optimizations.compute import (submit_molecules,
                                             export_molecule_data,
                                             execute_optimization_from_server,
                                             execute_optimization_from_molecule)
                    


@click.group()
def cli():
    pass

@cli.group()
def fractal():
    pass

@fractal.command()
def server_init():
    print("Fractal init!")

@fractal.command()
def server_start():
    print("Fractal start!")

@fractal.command()
def manager_start():
    print("Fractal start!")

@cli.group()
def optimize():
    pass

@optimize.command()
def submit():
    pass

@optimize.command()
def export():
    pass

@optimize.command()
def status():
    pass

@optimize.command()
@click.option('--sleep-time')
def error_cycle():
    pass

@optimize.command()
def local_from_server():
    pass

@optimize.command()
def local_from_molecule():
    pass
