"""
This module defines the layout and style of the webapp
"""

import dash_html_components as html
import dash_core_components as dcc
import dash_bootstrap_components as dbc
import dash_bio

from .content import available_datasets, available_forcefields, FORCEFIELDS
from .plots import PLOTS

THEME = dbc.themes.JOURNAL


def _header():
    """
    This function creates the top part of the webapp
    """
    return (
        dbc.Alert(
            "The app is in an alpha state and may contain incorrect results.", color="warning"
        ),
        html.H1("Benchmark Dataset Viewer"),
        html.Div("Useful benchmarks for the Openforcefield Initiative"),
        html.Hr(),
    )


def _controls():
    """
    This function creates all the widgets used to configure the plots.

    Components
    ----------
    * `controls-datasets`: Dataset dropdown
    * `controls-forcefields`: Forcefield dropdown
    * `controls-plots`: Which plot to graph
    * `controls-plot-config`: Configure the plot
    """
    return (
        dbc.Row(
            [
                dbc.Col([dbc.Label("Choose a dataset:")], width=2),
                dbc.Col(
                    [
                        dcc.Dropdown(
                            id="controls-datasets",
                            options=available_datasets(),
                            value=available_datasets()[0]["value"],
                        )
                    ]
                ),
                dbc.Col([dbc.Label("Choose forcefields:")], width=2),
                dbc.Col(
                    [
                        dcc.Dropdown(
                            id="controls-forcefields",
                            options=available_forcefields(),
                            value=list(FORCEFIELDS.keys()),
                            multi=True,
                        )
                    ]
                ),
            ]
        ),
        dbc.Row(
            [
                dbc.Col([dbc.Label("Choose a plot:")], width=2),
                dbc.Col(
                    [dcc.Dropdown(id="controls-plots", options=PLOTS, value=PLOTS[0]["value"])]
                ),
                dbc.Col([dbc.Label("Plot parameters:")], width=2),
                dbc.Col([dbc.Label("Nothing to show yet.")], id="controls-plot-config"),
            ]
        ),
    )


def _main_plot():
    """
    This function creates the main plot area.

    Components
    ----------
    * `primary-graph`: the plotly graph itself
    * `info-dataset-name`: bar header of the plot
    * `graph-toast-error`: display errors here
    """
    return (
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
                    id="loading-graph", children=[dcc.Graph(id="primary-graph")], type="default"
                ),
            ]
        ),
    )


def _molecule_viewer():
    """
    This function creates the main plot area.

    Components
    ----------
    * `available-molecules`: which molecules can be displayed, listed in the dropdown
    * `molecule-toast-error`: display errors here
    * `dash-bio-3d`: the molecule viewer
    """
    return (
        dbc.Card(
            [
                dbc.CardHeader("Molecule Explorer"),
                dcc.Dropdown(id="available-molecules", options=[], multi=False, className="p-2"),
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
                            id="dash-bio-3d", styles={}, modelData={"atoms": [], "bonds": []}
                        )
                    ],
                    type="default",
                ),
            ]
        ),
    )


def _info():
    """
    This function creates the main plot area.

    Components
    ----------
    * `dataset-information`: display more information on the loaded dataset
    """
    return (
        dbc.Card(
            [
                dbc.CardHeader("Dataset Information"),
                dbc.ListGroup(id="dataset-information", flush=True),
            ]
        ),
    )


def layout():
    """
    Concatenates all component groups as defined in the functions above. This is
    the final result that will be rendered in the webapp.
    """
    return dbc.Container([*_header(), *_controls(), *_main_plot()])
