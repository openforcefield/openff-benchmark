"""
This module defines the layout and style of the webapp
"""
from functools import partial
import dash_html_components as html
import dash_core_components as dcc
import dash_bootstrap_components as dbc
import dash_bio

from .content import available_datasets, available_forcefields, FORCEFIELDS
from .plots import boxplots_csv

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
    * # DISABLED - `controls-plots`: Which plot to graph
    * # DISABLED - `controls-plot-config`: Configure the plot
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
            ]
        ),
        dbc.Row(
            [
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
    )


def _general_errors():
    """
    Toast with an error message

    Components
    ----------
    * `toast-error`
    * `toast-error-message`
    """
    return (
        dbc.Toast(
            [html.P(id="toast-error-message")],
            id="toast-error",
            header="An error occured!",
            icon="danger",
            dismissable=True,
            is_open=False,
            style={"max-width": "100%"},
        ),
    )


def _quantum_chemistry():
    """
    This function creates the plot area for quantum chemistry comparisons:

    - Structural accuracy
    - Conformer energies

    Components
    ----------
    * `structural`: plotly graph for RMSD accuracy (+ their errors; check `_plot_area`)
    * `energies`: plotly graph for conformer energies (+ their errors; check `_plot_area`)
    * `quantum-chemistry-toggle`: button controlling open/closed state of plot area
    * `quantum-chemistry-collapse`: collapsible area that contains plots
    """
    return (
        dbc.Card(
            [
                dbc.CardHeader(
                    html.H2(
                        dbc.Button(
                            f"Quantum Chemistry", color="link", id=f"quantum-chemistry-toggle",
                        )
                    )
                ),
                dbc.Collapse(
                    [*_plot_area("structural"), *_plot_area("energies"),],
                    id="quantum-chemistry-collapse",
                ),
            ]
        ),
    )


def _physical_properties():
    """
    This function creates the plot area for comparisons of physical properties:

    - Liquid density (???)
    - Liquid density (mixture)

    Components
    ----------
    * `liquid`: plotly graph for RMSD accuracy (+ their errors; check `_plot_area`)
    * `liquid-mixture`: plotly graph for conformer energies (+ their errors; check `_plot_area`)
    * `physical-properties-toggle`: button controlling open/closed state of plot area
    * `physical-properties-collapse`: collapsible area that contains plots
    """
    return (
        dbc.Card(
            [
                dbc.CardHeader(
                    html.H2(
                        dbc.Button(
                            "Physical Properties", color="link", id="physical-properties-toggle",
                        )
                    )
                ),
                dbc.Collapse(
                    [
                        dbc.Toast(
                            [html.P("These two plots are for demo purposes only - no real data")],
                            header="Warning!",
                            icon="danger",
                            dismissable=True,
                            is_open=True,
                            style={"max-width": "100%"},
                        ),
                        *_plot_area("liquid-density"),
                        *_plot_area("liquid-density-mixture"),
                    ],
                    id="physical-properties-collapse",
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


def _plot_area(identifier):
    """
    Loading area for plot and some errors if needed

    Components
    ----------
    * `graph-{identifier}-error`: Box containing error info, if needed
    * `graph-{identifier}-error-message`: Error message
    * `graph-{identifier}-loading`: Loading layer
    * `{identifier}`: The target dcc.Graph component
    """
    return (
        dbc.Toast(
            [html.P(id=f"{identifier}-error-message")],
            id=f"{identifier}-error",
            header="An error occured!",
            icon="danger",
            dismissable=True,
            is_open=False,
            style={"max-width": "100%"},
        ),
        dcc.Loading(
            id=f"{identifier}-loading", children=[dcc.Graph(id=identifier)], type="default",
        ),
    )


def layout():
    """
    Concatenates all component groups as defined in the functions above. This is
    the final result that will be rendered in the webapp.
    """
    return dbc.Container(
        [
            *_header(),
            *_controls(),
            *_general_errors(),
            *_quantum_chemistry(),
            *_physical_properties(),
        ]
    )


PLOTS = {
    "liquid-density": {
        "title": "Liquid density",
        "yaxis": "Error (%)",
        "plotter": partial(boxplots_csv, selector=lambda df: df.RMSD.values),
    },
    "liquid-density-mixture": {
        "title": "Liquid density (mixture)",
        "yaxis": "Error (%)",
        "plotter": partial(boxplots_csv, selector=lambda df: df.qm_energy.values),
    },
    "structural": {
        "title": "Structural Accuracy (RMSD)",
        "yaxis": "Distance (Å)",
        "plotter": partial(boxplots_csv, selector=lambda df: df.TFD.values),
    },
    "energies": {
        "title": "Conformer energies",
        "yaxis": "Energy (kcal/mol)",
        "plotter": partial(boxplots_csv, selector=lambda df: df.ff_energy.values),
    },
}
