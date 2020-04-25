import logging
import traceback

import dash
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate

from .content import get_dataframe_from_trim3, get_json_from_off, FORCEFIELDS
from .layout import PLOTS


logger = logging.getLogger(__name__)


def register_callbacks(app):

    figure_outputs = [
        [
            Output(figure_id, "figure"),
            Output(f"{figure_id}-error", "is_open"),
            Output(f"{figure_id}-error-message", "children"),
        ]
        for figure_id in sorted(PLOTS.keys())
    ]

    @app.callback(
        [
            Output("toast-error", "is_open"),
            Output("toast-error-message", "children"),
            *[out for figure in figure_outputs for out in figure],
        ],
        [Input("controls-datasets", "value"), Input("controls-forcefields", "value")],
    )
    def draw_plots(dataset, forcefields):
        """
        Creates the needed plots using functions in openbenchmark.plots.

        Returns
        -------
        tuple of length (2 + 3*n_figures)
            bool: True if error message should be displayed, False otherwise
            str or None: Error message, if any

            For each figure, append:
                dict: figure data for plotly Graph
                bool: True if error message should be displayed for this figure, False otherwise
                str or None: Error message for this figure, if any

        """
        if None in (dataset, forcefields):
            # Return is open/closed, error_msg
            return [False, None] + [{}, False, None] * len(PLOTS)
        # First, get all data
        try:
            jsondicts = list(get_dataframe_from_trim3(dataset, tuple(forcefields)))
            ff_names = [FORCEFIELDS[ff] for ff in forcefields]
        except Exception as exc:  # if error occured, print in general error box
            logger.error(traceback.format_exc())
            error_msg = "\n".join(traceback.format_exc(limit=0, chain=False).splitlines()[1:])
            return [True, error_msg] + [{}, False, None] * len(PLOTS)

        # Now, try to plot figures
        result = [False, None]
        for figure_id, meta in sorted(PLOTS.items()):
            try:
                figure = meta["plotter"](
                    jsondicts, forcefield_names=ff_names, title=meta["title"], yaxis=meta["yaxis"],
                )
                error = False
                error_msg = None
            except Exception as exc:
                logger.error(traceback.format_exc())
                figure = {}
                error = True
                error_msg = "\n".join(traceback.format_exc(limit=0, chain=False).splitlines()[1:])
            finally:
                result.extend([figure, error, error_msg])
        return result

    @app.callback(
        [
            Output(f"quantum-chemistry-collapse", "is_open"),
            Output(f"physical-properties-collapse", "is_open"),
        ],
        [
            Input(f"quantum-chemistry-toggle", "n_clicks"),
            Input(f"physical-properties-toggle", "n_clicks"),
        ],
        [
            State(f"quantum-chemistry-collapse", "is_open"),
            State(f"physical-properties-collapse", "is_open"),
        ],
    )
    def toggle_accordion(quantum, physical, quantum_isopen, physical_isopen):
        ctx = dash.callback_context

        if not ctx.triggered:
            return ""
        else:
            button_id = ctx.triggered[0]["prop_id"].split(".")[0]
        # Return state (open=, closed=False) for:
        # 1) Quantum Chemistry plots
        # 2) Physical properties plots
        if button_id == "quantum-chemistry-toggle" and quantum:
            return not quantum_isopen, False
        elif button_id == "physical-properties-toggle" and physical:
            return False, not physical_isopen
        return True, True

    return app
