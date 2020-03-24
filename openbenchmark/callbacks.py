import logging
import traceback

from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate

from .content import get_dataframe, FORCEFIELDS
from . import plots


logger = logging.getLogger(__name__)


def register_callbacks(app):
    @app.callback(
        [
            Output("primary-graph", "figure"),
            Output("graph-toast-error", "is_open"),
            Output("graph-toast-error-message", "children"),
        ],
        [
            Input("controls-datasets", "value"),
            Input("controls-forcefields", "value"),
            Input("controls-plots", "value"),
            Input("controls-plot-config", "children"),
        ],
    )
    def draw_plot(dataset, forcefields, plot, config):
        if None in (dataset, forcefields, plot):
            return {}, False, None
        try:
            dataframes = get_dataframe(dataset, forcefields)
            ff_names = [FORCEFIELDS[ff] for ff in forcefields]
            fig = getattr(plots, plot)(dataframes, ff_names)
            return fig, False, None
        except Exception as exc:
            logger.error(traceback.format_exc())
            error_msg = "\n".join(traceback.format_exc(limit=0, chain=False).splitlines()[1:])
            return {}, True, error_msg

    return app
