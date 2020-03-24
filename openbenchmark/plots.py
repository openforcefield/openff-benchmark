def dde_vs_rmsd(dataframes, forcefields):
    for df in dataframes:
        print(df.head())
    return {
        "data": [
            {"x": [1, 2, 3], "y": [4, 1, 2], "type": "bar", "name": "SF"},
            {"x": [1, 2, 3], "y": [2, 4, 5], "type": "bar", "name": u"Montréal"},
        ],
        "layout": {"title": "Dash Data Visualization"},
    }


def _histogram(dataframes, forcefields, column, title):
    plot = {
        "data": [],
        "layout": {
            "title": title,
            "barmode": "overlay",
            "hovermode": "closest",
            "yaxis": {"title": {"text": "Frequency (%)"}},
        },
    }
    for df, ff in zip(dataframes, forcefields):
        plot["data"].append(
            {
                "x": df[column].values,
                "type": "histogram",
                "nbinsx": 250,
                "name": ff,
                "legendgroup": ff,
                "histnorm": "percent",
                "opacity": 0.5,
            }
        )
    return plot


def rmsd_distribution(dataframes, forcefields):
    return _histogram(dataframes, forcefields, column="RMSD", title="RMSD (Å)")


def tfd_distribution(dataframes, forcefields):
    return _histogram(dataframes, forcefields, column="TFD", title="Torsion Fingerprint Deviation")


def energy_mse(dataframes, forcefields):

    return {
        "data": [
            {"x": [1, 2, 3], "y": [4, 1, 2], "type": "bar", "name": "SF"},
            {"x": [1, 2, 3], "y": [2, 4, 5], "type": "bar", "name": u"Montréal"},
        ],
        "layout": {"title": "Dash Data Visualization"},
    }


PLOTS = (
    {"label": "RMSD distribution", "value": "rmsd_distribution"},
    {"label": "TFD distribution", "value": "tfd_distribution"},
    # {"label": "ddE vs RMSD", "value": "dde_vs_rmsd"},
    # {"label": "Energy Mean Signed Errors", "value": "energy_mse"},
)
