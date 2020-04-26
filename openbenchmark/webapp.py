"""
webapp.py
Comparison benchmarks between public force fields and Open Force Field Initiative force fields
"""

from dash import Dash

from .layout import layout, THEME
from .callbacks import register_callbacks


def create_app():
    """
    Bootstrap Dash app with layout, theme and callbacks
    """

    app = Dash(__name__, external_stylesheets=[THEME])
    app.title = "Openbenchmark - Openforcefield Initiative"
    app.layout = layout()
    app = register_callbacks(app)
    return app


def main():
    app = create_app()
    app.run_server(debug=True)


if __name__ == "__main__":
    main()
