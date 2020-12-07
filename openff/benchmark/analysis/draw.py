#! /usr/bin/env python

"""
draw.py

Plot generation for openff-benchmark

By:      David F. Hahn, Victoria T. Lim
Version: Dec 7 2020

"""

import os
import numpy as np
import pandas as pd
from scipy.interpolate import interpn
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from plotly import graph_objects as go

def draw_scatter(
    x_data, y_data, method_label, x_label, y_label, out_file, what_for="talk"
):
    """
    Draw scatter plot, such as of (ddE vs RMSD) or (ddE vs TFD).

    Parameters
    ----------
    x_data : list
        x_data[i] represents ith molecular structure
    y_data : list
        should have same shape and correspond to x_data
    method_label : str
        list of all the method names including reference method first
    x_label : string
        name of the x-axis label
    y_label : string
        name of the y-axis label
    out_file : string
        name of the output file
    what_for : string
        dictates figure size, text size of axis labels, legend, etc.
        "paper" or "talk"

    """
    markers = ["o", "^", "d", "x", "s", "p", "P", "3", ">"]

    p = plt.scatter(
        x_data,
        y_data,
#        marker=markers[i],
        label=method_label,
        alpha=0.6,
    )

    if what_for == "paper":
        fig = plt.gcf()
        fig.set_size_inches(4, 3)
        plt.subplots_adjust(left=0.16, right=0.72, top=0.9, bottom=0.2)
        plt.xlabel(x_label, fontsize=10)
        plt.ylabel(y_label, fontsize=10)
        plt.xticks(fontsize=10)
        plt.yticks(fontsize=10)
        plt.legend(loc=(1.04, 0.4), fontsize=10)
        # make the marker size smaller
        p.set_sizes([8.0])

    elif what_for == "talk":
        plt.xlabel(x_label, fontsize=14)
        plt.ylabel(y_label, fontsize=14)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.legend(loc=(1.04, 0.5), fontsize=12)
        # make the marker size smaller
        p.set_sizes([4.0])

    # set log scaling but use symmetric log for negative values
    #    plt.yscale('symlog')

    plt.savefig(out_file, bbox_inches="tight")
    plt.clf()
    # plt.show()



def draw_corr(
    x_data, y_data, method, x_label, y_label, out_file, what_for="talk"
):
    """
    Draw scatter plot, such as of (ddE vs RMSD) or (ddE vs TFD).

    Parameters
    ----------
    x_data : list of lists
        x_data[i][j] represents ith method and jth molecular structure
    y_data : list of lists
        should have same shape and correspond to x_data
    method_labels : list
        list of all the method names including reference method first
    x_label : string
        name of the x-axis label
    y_label : string
        name of the y-axis label
    out_file : string
        name of the output file
    what_for : string
        dictates figure size, text size of axis labels, legend, etc.
        "paper" or "talk"

    """
    p = plt.scatter(
        x_data,
        y_data,
        label=method,
        alpha=0.6,
    )

    if what_for == "paper":
        fig = plt.gcf()
        fig.set_size_inches(5, 4)
        plt.subplots_adjust(left=0.16, right=0.72, top=0.9, bottom=0.2)
        plt.xlabel(x_label, fontsize=10)
        plt.ylabel(y_label, fontsize=10)
        plt.xticks(fontsize=10)
        plt.yticks(fontsize=10)
        plt.legend(loc=(1.04, 0.4), fontsize=10)
        # make the marker size smaller
        p.set_sizes([8.0])

    elif what_for == "talk":
        plt.xlabel(x_label, fontsize=14)
        plt.ylabel(y_label, fontsize=14)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.legend(loc=(1.04, 0.5), fontsize=12)
        # make the marker size smaller
        p.set_sizes([4.0])

    plt.savefig(out_file, bbox_inches="tight")
    plt.clf()
    # plt.show()


def draw_ridgeplot(
    dataframes,
    key,
    x_label,
    out_file,
    what_for="paper",
    bw="hist",
    same_subplot=False,
    sym_log=False,
    hist_range=(-15, 15)
):

    """
    Draw ridge plot of data (to which kernel density estimate is applied)
    segregated by each method (representing a different color/level).

    Modified from the following code:
    https://seaborn.pydata.org/examples/kde_ridgeplot.html

    Parameters
    ----------
    dataframe : dict of pandas.DataFrame 
    key : string
        key specifying the variable to be plotted
    x_label : string
        name of the x-axis label
        also used for pandas dataframe column name
    out_file : string
        name of the output file
    what_for : string
        dictates figure size, text size of axis labels, legend, etc.
        "paper" or "talk"
    bw : string or float
        defines bandwidth for KDE as called in seaborn.kdeplot, OR don't use
        kde at all and just histogram the data;
        options: 'scott' (KDE), 'silverman' (KDE), scalar value (KDE), 'hist'
    same_subplot : Boolean
        False is default to have separate and slightly overlapping plots,
        True to plot all of them showing on the same subplot (no fill)
    sym_log : Boolean
        False is default to plot density estimate as is,
        True to plot x-axis on symmetric log scale
    hist_range : tuple
        tuple of min and max values to use for histogram;
        only needed if bw is set to 'hist'

    """
    if what_for == "paper":
        ridgedict = {
            "h": 0.9,
            "lw": 2.0,
            "vl": 1.0,
            "xfontsize": 14,
        }
    elif what_for == "talk":
        ridgedict = {
            "h": 2.0,
            "lw": 3.0,
            "vl": 1.0,
            "xfontsize": 16,
        }

    # Initialize the FacetGrid object
    my_cmap = "tab10"
    sns.palplot(sns.color_palette(my_cmap))
    colors = sns.color_palette(my_cmap)


    fig = go.Figure()
    for method, result in dataframes.items():
        if bw == "hist":
            hist = np.histogram(result[key], bins=15, range=hist_range, density=True)
            fig.add_trace(
                go.Scatter(
                    x=hist[1],
                    y=hist[0]+[hist[0][-1]],
                    name=method,
                    mode='lines',
                    line=dict(
                        width=ridgedict["lw"],
                        shape='hv'
                    ),
                )
            )
        elif bw == "kde":
            fig.add_trace(
                go.Violin(
                    x=result[key],
                    line_width=ridgedict["lw"]
                )
            )
            

    if bw == "kde":
        fig.update_traces(
            orientation="h",
            side="positive",
            points=False
            )

    fig.update_layout(
        xaxis=dict(
            title=dict(
                text=x_label,
                font_size=ridgedict["xfontsize"]
                ),
            range=hist_range
            ),
        shapes=[
            dict(
                type= 'line',
                yref= 'paper', y0= 0, y1= 1,
                xref= 'x', x0= 0, x1= 0
            )
        ]
    )
    
    fig.write_image(out_file)
    

def draw_density2d(
    x_data,
    y_data,
    title,
    x_label,
    y_label,
    out_file,
    what_for="talk",
    bins=20,
    x_range=None,
    y_range=None,
    z_range=None,
    z_interp=True,
    symlog=False,
):

    """
    Draw a scatter plot colored smoothly to represent the 2D density.
    Based on: https://stackoverflow.com/a/53865762/8397754

    Parameters
    ----------
    x_data : 1D list
        represents x-axis data for all molecules of a given method
    y_data : 1D list
        should have same shape and correspond to x_data
    title : string
        title of the plot
    x_label : string
        name of the x-axis label
    y_label : string
        name of the y-axis label
    out_file : string
        name of the output file
    what_for : string
        dictates figure size, text size of axis labels, legend, etc.
        "paper" or "talk"
    bins : int
        number of bins for np.histogram2d
    x_range : tuple of two floats
        min and max values of x-axis
    y_range : tuple of two floats
        min and max values of y-axis
    z_range : tuple of two floats
        min and max values of density for setting a uniform color bar;
        these should be at or beyond the bounds of the min and max
    z_interp : Boolean
        True to smoothen the color scale for the scatter plot points;
        False to plot 2d histograms colored by cells (no scatter plot)
    symlog : Boolean
        True to represent y-axis on (symmetric) log scale (linear
        between -1 and 1 ), False for linear y-scaling
    """

    def colorbar_and_finish(labelsize, fname):
        cb = plt.colorbar()
        cb.ax.tick_params(labelsize=labelsize)
        cb.ax.set_title("counts", size=labelsize)

        plt.savefig(fname, bbox_inches="tight")
        plt.clf()
        # plt.show()

    fig = plt.gcf()
    if what_for == "paper":
        ms = 1
        size1 = 10
        size2 = 10
        fig.set_size_inches(4, 3)
    elif what_for == "talk":
        ms = 4
        size1 = 14
        size2 = 16
        fig.set_size_inches(9, 6)
    plt_options = {"s": ms, "cmap": "coolwarm_r"}

    # label and adjust plot
    plt.title(title, fontsize=size2)
    plt.xlabel(x_label, fontsize=size2)
    plt.ylabel(y_label, fontsize=size2)
    plt.xticks(fontsize=size1)
    plt.yticks(fontsize=size1)

    if x_range is not None:
        plt.xlim(x_range[0], x_range[1])

    if y_range is not None:
        plt.ylim(y_range[0], y_range[1])

    # remove any nans from x_data, such as TFD score for urea-like mols
    nan_inds = np.isnan(x_data)
    x_data = x_data[~nan_inds]
    y_data = y_data[~nan_inds]
    print(f"\nNumber of data points in FF scatter plot: {len(x_data)}")

    # compute histogram in 2d
    data, x_e, y_e = np.histogram2d(x_data, y_data, bins=bins)

    # plot colored 2d histogram if z_interp not specified
    if not z_interp:
        extent = [x_e[0], x_e[-1], y_e[0], y_e[-1]]
        plt.imshow(
            data.T,
            extent=extent,
            origin="lower",
            aspect="auto",
            cmap=plt_options["cmap"],
            vmin=z_range[0],
            vmax=z_range[1],
        )

        colorbar_and_finish(size1, out_file)
        return

    # smooth/interpolate data
    z = interpn(
        (0.5 * (x_e[1:] + x_e[:-1]), 0.5 * (y_e[1:] + y_e[:-1])),
        data,
        np.vstack([x_data, y_data]).T,
        method="splinef2d",
        bounds_error=False,
    )

    # sort the points by density, so that the densest points are plotted last
    idx = z.argsort()
    x, y, z = x_data[idx], y_data[idx], z[idx]

    print(
        f"{title} ranges of data in density plot:\n\t\tmin\t\tmax"
        f"\nx\t{np.min(x):10.4f}\t{np.max(x):10.4f}"
        f"\ny\t{np.min(y):10.4f}\t{np.max(y):10.4f}"
        f"\nz\t{np.min(data):10.4f}\t{np.max(data):10.4f} (histogrammed)"
        f"\nz'\t{np.nanmin(z):10.4f}\t{np.nanmax(z):10.4f} (interpolated)"
    )

    # add dummy points of user-specified min/max z for uniform color scaling
    # similar to using vmin/vmax in pyplot pcolor
    x = np.append(x, [-1, -1])
    y = np.append(y, [0, 0])
    z = np.append(z, [z_range[0], z_range[1]])

    print(f"z''\t{np.nanmin(z):10.4f}\t{np.nanmax(z):10.4f} (interp, bounded)")

    # generate the plot
    plt.scatter(x, y, c=z, vmin=z_range[0], vmax=z_range[1], **plt_options)

    # set log scaling but use symmetric log for negative values
    if symlog:
        plt.yscale("symlog")

    # configure color bar and finish plotting
    colorbar_and_finish(size1, out_file)


def plot_compare_ffs(results_dir):
    results = {}
    for path in results_dir:
        if (os.path.isfile(path) and path.split('.')[-1].lower() == 'csv'):
            method = '.'.join(os.path.split(path)[1].split('.')[:-1])
            results[method] = pd.read_csv(path)
        elif os.path.isdir(path):
            for root, dirs, files in os.walk(path):
                for file in files:
                    path = os.path.join(root, file)
                    if (os.path.isfile(path) and path.split('.')[-1].lower() == 'csv'):
                        method = '.'.join(file.split('.')[:-1])
                        results[method] = pd.read_csv(path)

    for method, result in results.items():
        # ddE vs RMSD scatter
        draw_scatter(result.loc[:, 'rmsd'],
                     result.loc[:,'dde'],
                     method,
                     'rmsd',
                     'dde',
                     f'fig_{method}_scatter_rmsd_dde.png',
                     what_for="talk"
                     )
        draw_scatter(result.loc[:, 'tfd'],
                     result.loc[:,'dde'],
                     method,
                     'tfd',
                     'dde',
                     f'fig_{method}_scatter_tfd_dde.png',
                     what_for="talk"
                     )
        
        # draw density 2D
        draw_density2d(
            result.loc[:, 'rmsd'],
            result.loc[:,'dde'],
            method,
            'rmsd',
            'dde',
            f'fig_{method}_density_rmsd_dde_linear.png',
            what_for="talk",
            x_range=(0,3.7),
            y_range=(-30,30),
            z_range=(-260,5200),
            z_interp=True,
            symlog=False
        )
        draw_density2d(
            result.loc[:, 'tfd'],
            result.loc[:,'dde'],
            method,
            'tfd',
            'dde',
            f'fig_{method}_density_rmsd_dde_linear.png',
            what_for="talk",
            x_range=(0,0.8),
            y_range=(-30,30),
            z_range=(-302,7060),
            z_interp=True,
            symlog=False
        )

        draw_density2d(
            result.loc[:,'rmsd'],
            result.loc[:,'dde'],
            method,
            'rmsd',
            'dde',
            f'fig_{method}_density_rmsd_dde_log.png',
            what_for="talk",
            x_range=(0,3.7),
            y_range=(-30,30),
            z_range=(-260,5200),
            z_interp=True,
            symlog=True
        )
        draw_density2d(
            result.loc[:, 'tfd'],
            result.loc[:,'dde'],
            method,
            'tfd',
            'dde',
            f'fig_{method}_density_rmsd_dde_log.png',
            what_for="talk",
            x_range=(0,0.8),
            y_range=(-30,30),
            z_range=(-302,7060),
            z_interp=True,
            symlog=True
        )

        # draw correlation plots between methods
        for method2 in list(results.keys())[list(results.keys()).index(method)+1:]:
            draw_corr(
                result.loc[:, 'rmsd'],
                results[method2].loc[:, 'rmsd'],
                f"{method2} vs. {method}",
                f"RMSD {method} ($\AA$)",
                f"RMSD {method2} ($\AA$)",
                f"fig_scatter_rmsd_{method2}_{method}.png",
                "paper"
            )

            draw_corr(
                result.loc[:, 'tfd'],
                results[method2].loc[:, 'tfd'],
                f"{method2} vs. {method}",
                f"TFD {method}",
                f"TFD {method2}",
                f"fig_scatter_tfd_{method2}_{method}.png",
                "paper"
            )

            draw_corr(
                result.loc[:, 'dde'],
                results[method2].loc[:, 'dde'],
                f"{method2} vs. {method}",
                f"ddE {method} (kcal/mol)",
                f"ddE {method2} (kcal/mol)",
                f"fig_scatter_dde_{method2}_{method}.png",
                "paper"
            )
            
        
    draw_ridgeplot(
        results,
        'dde',
        'ddE (kcal/mol)',
        out_file=f'fig_ridge_dde.png',
        what_for="talk",
        bw="kde",
        same_subplot=True,
        sym_log=False,
        hist_range=(-15,15)
    )
    draw_ridgeplot(
        results,
        'rmsd',
        'RMSD ($\AA$)',
        out_file=f'fig_ridge_rmsd.png',
        what_for="talk",
        bw="hist",
        same_subplot=True,
        sym_log=False,
        hist_range=(0, 3)
    )
    draw_ridgeplot(
        results,
        'tfd',
        'TFD',
        out_file=f'fig_ridge_tfd.png',
        what_for="talk",
        bw="hist",
        same_subplot=True,
        sym_log=False,
        hist_range=(0, 0.5)
    )
    



