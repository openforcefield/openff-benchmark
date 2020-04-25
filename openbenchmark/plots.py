import numpy as np
from functools import partial


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


def boxplots(dicts, forcefield_names=None, title="Boxplot", yaxis="Error"):
    all_data = []
    for d in dicts:
        data = {"y": np.random.rand(50), "type": "box", "name": d["forcefield"]}
        all_data.append(data)
    return {
        "data": all_data,
        "layout": {"title": title, "yaxis": {"title": {"text": yaxis}},},
    }


def boxplots_csv(dfs, selector, forcefield_names=None, title="Boxplot", yaxis="Error"):
    all_data = []
    for df, ff in zip(dfs, forcefield_names):
        data = {"y": selector(df), "type": "box", "name": ff}
        all_data.append(data)
    return {
        "data": all_data,
        "layout": {"title": title, "yaxis": {"title": {"text": yaxis}},},
    }
