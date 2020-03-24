"""
This module defines the layout and style of the webapp
"""

import dash_html_components as html
import dash_core_components as dcc
import dash_bootstrap_components as dbc
import dash_bio

THEME = dbc.themes.JOURNAL


def layout():
    return dbc.Container(
        [
            # Header
            dbc.Alert(
                "The app is in an alpha state and may contain incorrect results.", color="warning"
            ),
            html.H1("Benchmark Dataset Viewer"),
            html.Div("Useful benchmarks for the Openforcefield initiative"),
            # Primary data visualizer
            dbc.Card(
                [
                    dbc.CardHeader(id="info-dataset-name"),
                    dbc.Toast(
                        [html.P(id="graph-toast-error-message")],
                        id="graph-toast-error",
                        header="An error occured!",
                        icon="danger",
                        dismissable=True,
                        is_open=False,
                        style={"max-width": "100%"},
                    ),
                    dcc.Loading(
                        id="loading-graph",
                        # Callbacks can target this dcc.Graph using id=primary-graph
                        children=[dcc.Graph(id="primary-graph")],
                        type="default",
                    ),
                    dbc.Label(children="Benchmark:", id="graph-benchmark-label"),
                ]
            ),
            ### Molecule Explorer
            dbc.Card(
                [
                    dbc.CardHeader("Molecule Explorer"),
                    dcc.Dropdown(
                        id="available-molecules", options=[], multi=False, className="p-2"
                    ),
                    dbc.Toast(
                        [html.P(id="molecule-toast-error-message")],
                        id="molecule-toast-error",
                        header="An error occured!",
                        icon="danger",
                        dismissable=True,
                        is_open=False,
                        style={"max-width": "100%"},
                    ),
                    dcc.Loading(
                        id="loading-molecule",
                        children=[
                            dash_bio.Molecule3dViewer(
                                # Callbacks can target this 3D viewer using id=dash-bio-3d
                                id="dash-bio-3d",
                                styles={},
                                modelData={"atoms": [], "bonds": []},
                            )
                        ],
                        type="default",
                    ),
                ]
            ),
            ### Info
            dbc.Card(
                [
                    dbc.CardHeader("Dataset Information"),
                    # Callbacks can target this panel using id=dataset-information
                    dbc.ListGroup(id="dataset-information", flush=True),
                ]
            ),
        ]
    )
