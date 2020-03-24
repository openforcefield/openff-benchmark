"""
webapp.py
Comparison benchmarks between public force fields and Open Force Field Initiative force fields
"""


def create_app():
    """
    Bootstrap Dash app with layout, theme and callbacks
    """
    import dash
    from .layout import layout, THEME

    app = dash.Dash(__name__, external_stylesheets=[THEME])
    app.layout = layout()
    return app


def main():
    app = create_app()
    app.run_server(debug=True)


if __name__ == "__main__":
    main()
