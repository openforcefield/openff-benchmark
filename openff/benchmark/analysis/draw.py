#! /usr/bin/env python

"""
draw.py

Plot generation for openff-benchmark

By:      David F. Hahn, Lorenzo D'Amore, Victoria T. Lim
Version: Jul 22 2021

"""

import os
import numpy as np
import pandas as pd
from scipy.interpolate import interpn
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import warnings


def get_results(results_dir):
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

    return results, method


def res_inters(results, method, idx):
    """
    Takes only conformer which are present in all methods.

    Modifies `results` dictionary Dataframes in-place.

    Parameters
    ----------
    results : Dict[str, pd.DataFrame]
        Dict of Dataframes with method as keys, results as values.
        This function modifies the Dataframes in-place.
    method : str
        Key in results used as the first method for intersection.
        Can ultimately be any key in the dictionary due to the nature of the intersection.
    idx : string
        idx is the column of the results df containing the molecule index:
        'name' for compare-ff, match-minima
        'mm_ref' for lucas analysis
        'mm_conf' for swope analysis
    """
    for m, df in results.items():
        results[m].set_index(idx, inplace=True)

    index_intersect = results[method].index
    for m in results:
        index_intersect = index_intersect.intersection(results[m].index)
    for m, df in results.items():
        results[m] = df.loc[index_intersect]
        if results[m].shape != df.shape:
            warnings.warn(f"Not all conformers of method {m} considered, because these are not available in other methods.")


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
    y_label,
    out_file,
    what_for="paper",
    bw="hist",
    same_subplot=False,
    sym_log=False,
    hist_range=(-15, 15),
    stat="probability"
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
    stat : string
        normalization for the distribution, either "probability" or "count"

    """
    # Define and use a simple function to label the plot in axes coordinates
    def label(x, color, label):
        ax = plt.gca()

        # set axis font size
        if what_for == "paper":
            fs = 14
        elif what_for == "talk":
            fs = 14

        ax.text(
            0,
            0.2,
            label,
            fontweight="bold",
            color=color,
            fontsize=fs,
            ha="left",
            va="center",
            transform=ax.transAxes,
        )

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

    num_methods = len(dataframes)

    # Initialize the FacetGrid object
    my_cmap = "tab10"
    sns.palplot(sns.color_palette(my_cmap))
    colors = sns.color_palette(my_cmap)

    temp = []
    for method, result in dataframes.items():
        result = result.rename(columns={key: x_label})
        result["method"] = method
        temp.append(result)

    # list of dataframes concatenated to single dataframe
    df = pd.concat(temp, ignore_index=True)
    g = sns.FacetGrid(
        df, row="method", hue="method", aspect=10, height=ridgedict["h"], palette=colors
    )

    if not same_subplot:

        # draw filled-in densities
        if bw == "hist":
            histoptions = {
                "histtype": "bar",
                "alpha": 0.6,
                "linewidth": ridgedict["lw"],
                "range": hist_range,
                "align": "mid",
            }
            g.map(
                sns.displot,
                x_label,
                kind="hist",
                bins=15,
                stat=stat,
                hist_kws=histoptions,
            )
        else:
            g.map(
                sns.kdeplot,
                x_label,
                clip_on=False,
                shade=True,
                alpha=0.5,
                lw=ridgedict["lw"],
                bw=bw,
            )

        # draw colored horizontal line below densities
        g.map(plt.axhline, y=0, lw=ridgedict["lw"], clip_on=False)

    else:

        # draw black horizontal line below densities
        plt.axhline(y=0, color="black")

    # draw outline around densities; can also single outline color: color="k"
    if bw == "hist":
        histoptions = {
            "element": "step",
            "alpha": 1.0,
            "linewidth": ridgedict["lw"],
            "fill": False,
            "binrange": hist_range,
        }
        g.map(
            sns.histplot,
            x_label,
            bins=15,
            stat=stat,
            **histoptions,
        )

    else:
        g.map(sns.kdeplot, x_label, clip_on=False, lw=ridgedict["lw"], bw=bw)

    # draw a vertical line at x=0 for visual reference
    g.map(plt.axvline, x=0, lw=ridgedict["vl"], ls="--", color="gray", clip_on=False)

    # optional: add custom vertical line
    # g.map(plt.axvline, x=0.12, lw=1, ls='--', color='gray', clip_on=False)

    # add labels to each level
    if not same_subplot:
        g.map(label, x_label)

    # else if single subplot, generate a custom legend
    else:
        cmap = mpl.cm.tab10
        patches = []
        for i, label in enumerate(dataframes):
            patches.append(
                mpl.patches.Patch(
                    color=cmap(i/10),
                    label=label,
                )
            )
        plt.legend(handles=patches, fontsize=ridgedict["xfontsize"] / 1.2)

    # optional: set symmetric log scale on x-axis
    if sym_log:
        g.set(xscale="symlog")

    # Set the subplots to overlap
    if not same_subplot:
        g.fig.subplots_adjust(hspace=-0.45)
    else:
        g.fig.subplots_adjust(hspace=-1.0)

    # Remove axes details that don't play well with overlap
    g.set_titles("")
    #    g.set(yticks=[])
    g.despine(bottom=True)  # , left=True)
    # ax = plt.gca()
    # ax.spines['left'].set_visible(True)
    # ax.spines['left'].set_position('zero')
    # ax.set_yticks([0.4])
    if what_for == "paper":
        plt.gcf().set_size_inches(7, 3)
    elif what_for == "talk":
        plt.gcf().set_size_inches(12, 9)

    # adjust font sizes
    plt.xlabel(x_label, fontsize=ridgedict["xfontsize"])
    plt.ylabel(y_label, fontsize=ridgedict["xfontsize"])
    plt.xticks(fontsize=ridgedict["xfontsize"])

    # save with transparency for overlapping plots
    plt.savefig(out_file, transparent=True, bbox_inches="tight")
    plt.clf()


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

    if x_range is None:
        x_range = (0, x_data.max())
    plt.xlim(x_range[0], x_range[1])
    if y_range is None:
        yabs = y_data.abs().max()
        y_range = (-yabs, yabs)
    plt.ylim(y_range[0], y_range[1])

    # remove any nans from x_data, such as TFD score for urea-like mols
    nan_inds = x_data.isna()
    x_data = x_data.dropna().values
    y_data = y_data[~nan_inds].values

    # compute histogram in 2d
    data, x_e, y_e = np.histogram2d(x_data, y_data, bins=bins)
 #   print('data', data, data.min(), data.max())
 #   print('e', x_e, y_e)
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
#    print(z, z.max(), z.min())

    # sort the points by density, so that the densest points are plotted last
    idx = z.argsort()
    x, y, z = x_data[idx], y_data[idx], z[idx]


    if z_range is None:
        z_range=(0, np.max(z))
    # # add dummy points of user-specified min/max z for uniform color scaling
    # # similar to using vmin/vmax in pyplot pcolor
    x = np.append(x, [-1, -1])
    y = np.append(y, [0, 0])
    z = np.append(z, [z_range[0], z_range[1]])


    # generate the plot
    plt.scatter(x, y, c=z, vmin=z_range[0], vmax=z_range[1], **plt_options)

    # set log scaling but use symmetric log for negative values
    if symlog:
        plt.yscale("symlog")

    # configure color bar and finish plotting
    colorbar_and_finish(size1, out_file)


def plot_compare_ffs(results_dir, output_directory):
    os.makedirs(output_directory, exist_ok=True)

    results, method = get_results(results_dir)

    res_inters(results, method, 'name')

    plot_mol_minima(results, 'dde(kcal/mol)', 'dde[kcal/mol]', method, out_file=os.path.join(output_directory, 'minimaE.png'))
    
    plot_violin_signed(results, 'dde[kcal/mol]', out_file=os.path.join(output_directory, 'violin.svg'))


    for method, result in results.items():
        # ddE vs RMSD scatter
        draw_scatter(result.loc[:, 'rmsd'],
                     result.loc[:,'dde[kcal/mol]'],
                     method,
                     'rmsd',
                     'dde[kcal/mol]',
                     os.path.join(output_directory, f'fig_{method}_scatter_rmsd_dde.png'),
                     what_for="talk"
                     )
        draw_scatter(result.loc[:, 'tfd'],
                     result.loc[:,'dde[kcal/mol]'],
                     method,
                     'tfd',
                     'dde[kcal/mol]',
                     os.path.join(output_directory, f'fig_{method}_scatter_tfd_dde.png'),
                     what_for="talk"
                     )
        
        # draw density 2D
        draw_density2d(
            result.loc[:, 'rmsd'],
            result.loc[:,'dde[kcal/mol]'],
            method,
            'rmsd',
            'dde[kcal/mol]',
            os.path.join(output_directory, f'fig_{method}_density_rmsd_dde_linear.png'),
            what_for="talk",
            # x_range=(0,3.7),
            # y_range=(-30,30),
            # z_range=(-260,5200),
            z_interp=True,
            symlog=False
        )
        draw_density2d(
            result.loc[:, 'tfd'],
            result.loc[:,'dde[kcal/mol]'],
            method,
            'tfd',
            'dde[kcal/mol]',
            os.path.join(output_directory, f'fig_{method}_density_rmsd_dde_linear.png'),
            what_for="talk",
            # x_range=(0,0.8),
            # y_range=(-30,30),
            # z_range=(-302,7060),
            z_interp=True,
            symlog=False
        )

        draw_density2d(
            result.loc[:,'rmsd'],
            result.loc[:,'dde[kcal/mol]'],
            method,
            'rmsd',
            'dde[kcal/mol]',
            os.path.join(output_directory, f'fig_{method}_density_rmsd_dde_log.png'),
            what_for="talk",
            # x_range=(0,3.7),
            # y_range=(-30,30),
            # z_range=(-260,5200),
            z_interp=True,
            symlog=True
        )
        draw_density2d(
            result.loc[:, 'tfd'],
            result.loc[:,'dde[kcal/mol]'],
            method,
            'tfd',
            'dde[kcal/mol]',
            os.path.join(output_directory, f'fig_{method}_density_rmsd_dde_log.png'),
            what_for="talk",
            # x_range=(0,0.8),
            # y_range=(-30,30),
            # z_range=(-302,7060),
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
                os.path.join(output_directory, f"fig_scatter_rmsd_{method2}_{method}.png"),
                "paper"
            )

            draw_corr(
                result.loc[:, 'tfd'],
                results[method2].loc[:, 'tfd'],
                f"{method2} vs. {method}",
                f"TFD {method}",
                f"TFD {method2}",
                os.path.join(output_directory, f"fig_scatter_tfd_{method2}_{method}.png"),
                "paper"
            )

            draw_corr(
                result.loc[:, 'dde[kcal/mol]'],
                results[method2].loc[:, 'dde[kcal/mol]'],
                f"{method2} vs. {method}",
                f"ddE {method} (kcal/mol)",
                f"ddE {method2} (kcal/mol)",
                os.path.join(output_directory, f"fig_scatter_dde_{method2}_{method}.png"),
                "paper"
            )
            
        
    draw_ridgeplot(
        results,
        'dde[kcal/mol]',
        'ddE (kcal/mol)',
        'Density',
        out_file=os.path.join(output_directory, f'fig_ridge_dde.png'),
        what_for="talk",
        bw="hist",
        same_subplot=True,
        sym_log=False,
        hist_range=(-15,15),
        stat="probability"
    )
    draw_ridgeplot(
        results,
        'rmsd',
        'RMSD ($\AA$)',
        'Density',
        out_file=os.path.join(output_directory, f'fig_ridge_rmsd.png'),
        what_for="talk",
        bw="hist",
        same_subplot=True,
        sym_log=False,
        hist_range=(0, 3),
        stat="probability"
    )
    draw_ridgeplot(
        results,
        'tfd',
        'TFD',
        'Density',
        out_file=os.path.join(output_directory, f'fig_ridge_tfd.png'),
        what_for="talk",
        bw="hist",
        same_subplot=True,
        sym_log=False,
        hist_range=(0, 0.5),
        stat="probability"
    )
    

def plot_lucas(results_dir, output_directory):
    os.makedirs(output_directory, exist_ok=True)

    results, method = get_results(results_dir)

    res_inters(results, method, 'mm_ref')

    plot_violin_signed(results, 'dE (mm_ref-mm_min)', out_file=os.path.join(output_directory, 'violin.svg'))

    for method, result in results.items():
        # ddE vs RMSD scatter
        draw_scatter(result.loc[:, 'rmsd (mm_ref/mm_min)'],
                     result.loc[:,'dE (mm_ref-mm_min)'],
                     method,
                     'rmsd ($\AA$)',
                     'dE (kcal/mol)',
                     os.path.join(output_directory, f'fig_{method}_scatter_rmsd_de.png'),
                     what_for="talk"
                     )

        # draw density 2D
        draw_density2d(
            result.loc[:, 'rmsd (mm_ref/mm_min)'],
            result.loc[:,'dE (mm_ref-mm_min)'],
            method,
            'rmsd ($\AA$)',
            'dE (kcal/mol)',
            os.path.join(output_directory, f'fig_{method}_density_rmsd_de_linear.png'),
            what_for="talk",
            z_interp=True,
            symlog=False
        )

        draw_density2d(
            result.loc[:, 'rmsd (mm_ref/mm_min)'],
            result.loc[:,'dE (mm_ref-mm_min)'],
            method,
            'rmsd ($\AA$)',
            'dE (kcal/mol)',
            os.path.join(output_directory, f'fig_{method}_density_rmsd_de_log.png'),
            what_for="talk",
            z_interp=True,
            symlog=True
        )

        # draw correlation plots between methods
        for method2 in list(results.keys())[list(results.keys()).index(method)+1:]:
            draw_corr(
                result.loc[:, 'rmsd (mm_ref/mm_min)'],
                results[method2].loc[:, 'rmsd (mm_ref/mm_min)'],
                f"{method2} vs. {method}",
                f"rmsd {method} ($\AA$)",
                f"rmsd {method2} ($\AA$)",
                os.path.join(output_directory, f"fig_scatter_rmsd_{method2}_{method}.png"),
                "paper"
            )

            draw_corr(
                result.loc[:, 'dE (mm_ref-mm_min)'],
                results[method2].loc[:, 'dE (mm_ref-mm_min)'],
                f"{method2} vs. {method}",
                f"dE {method} (kcal/mol)",
                f"dE {method2} (kcal/mol)",
                os.path.join(output_directory, f"fig_scatter_de_{method2}_{method}.png"),
                "paper"
            )

    draw_ridgeplot(
        results,
        'dE (mm_ref-mm_min)',
        'dE (kcal/mol)',
        'Density',
        out_file=os.path.join(output_directory, f'fig_ridge_de.png'),
        what_for="talk",
        bw="hist",
        same_subplot=True,
        sym_log=False,
        hist_range=(-1.67,15),
        stat="probability"
    )

    draw_ridgeplot(
        results,
        'rmsd (mm_ref/mm_min)',
        'rmsd ($\AA$)',
        'Density',
        out_file=os.path.join(output_directory, f'fig_ridge_rmsd.png'),
        what_for="talk",
        bw="hist",
        same_subplot=True,
        sym_log=False,
        hist_range=(0, 3),
        stat="probability"
    )


def plot_swope(results_dir, de_cutoff, rmsd_cutoff, output_directory):
    os.makedirs(output_directory, exist_ok=True)

    results,method = get_results(results_dir)

    res_inters(results, method, 'mm_conf')

    plot_violin_signed(results, 'dE (mm_conf-mm_min)', out_file=os.path.join(output_directory, 'violin.svg'))

    for method, result in results.items():
        # ddE vs RMSD scatter
        draw_scatter(result.loc[:, 'rmsd (mm_conf/qm_min)'],
                     result.loc[:,'dE (mm_conf-mm_min)'],
                     method,
                     'rmsd ($\AA$)',
                     'dE (kcal/mol)',
                     os.path.join(output_directory, f'fig_{method}_scatter_rmsd_de.png'),
                     what_for="talk"
                     )

        # draw density 2D
        draw_density2d(
            result.loc[:, 'rmsd (mm_conf/qm_min)'],
            result.loc[:,'dE (mm_conf-mm_min)'],
            method,
            'rmsd ($\AA$)',
            'dE (kcal/mol)',
            os.path.join(output_directory, f'fig_{method}_density_rmsd_de_linear.png'),
            what_for="talk",
            z_interp=True,
            symlog=False
        )

        draw_density2d(
            result.loc[:, 'rmsd (mm_conf/qm_min)'],
            result.loc[:,'dE (mm_conf-mm_min)'],
            method,
            'rmsd ($\AA$)',
            'dE (kcal/mol)',
            os.path.join(output_directory, f'fig_{method}_density_rmsd_de_log.png'),
            what_for="talk",
            z_interp=True,
            symlog=True
        )

        # draw correlation plots between methods
        for method2 in list(results.keys())[list(results.keys()).index(method)+1:]:
            draw_corr(
                result.loc[:, 'rmsd (mm_conf/qm_min)'],
                results[method2].loc[:, 'rmsd (mm_conf/qm_min)'],
                f"{method2} vs. {method}",
                f"rmsd {method} ($\AA$)",
                f"rmsd {method2} ($\AA$)",
                os.path.join(output_directory, f"fig_scatter_rmsd_{method2}_{method}.png"),
                "paper"
            )

            draw_corr(
                result.loc[:, 'dE (mm_conf-mm_min)'],
                results[method2].loc[:, 'dE (mm_conf-mm_min)'],
                f"{method2} vs. {method}",
                f"dE {method} (kcal/mol)",
                f"dE {method2} (kcal/mol)",
                os.path.join(output_directory, f"fig_scatter_de_{method2}_{method}.png"),
                "paper"
            )

    combined_df = pd.concat(results, names=['method']).reset_index(0)

    new = combined_df[(combined_df['rmsd (mm_conf/qm_min)']<float(rmsd_cutoff))]
    new = new.reset_index()
    new.drop('mm_conf', axis='columns', inplace=True)
    new = new[['qm_min', 'mm_min', 'rmsd (mm_conf/qm_min)', 'dE (mm_conf-mm_min)', 'method']]
    new = { name : group for name,group in new.groupby('method')}

    new2 = combined_df[(combined_df['dE (mm_conf-mm_min)']<float(de_cutoff))]
    new2 = new2.reset_index()
    new2.drop('mm_conf', axis='columns', inplace=True)
    new2 = new2[['qm_min', 'mm_min', 'rmsd (mm_conf/qm_min)', 'dE (mm_conf-mm_min)', 'method']]
    new2 = { name : group for name,group in new2.groupby('method')}

    draw_ridgeplot(
        new,
        'dE (mm_conf-mm_min)',
        'dE (kcal/mol)',
        'Count',
        out_file=os.path.join(output_directory, f'fig_ridge_de_rmsd_cutoff_{rmsd_cutoff}.png'),
        what_for="talk",
        bw="hist",
        same_subplot=True,
        sym_log=False,
        hist_range=(-1.67,15),
        stat="count"
    )

    draw_ridgeplot(
        new2,
        'rmsd (mm_conf/qm_min)',
        'rmsd ($\AA$)',
        'Count',
        out_file=os.path.join(output_directory, f'fig_ridge_rmsd_de_cutoff_{de_cutoff}.png'),
        what_for="talk",
        bw="hist",
        same_subplot=True,
        sym_log=False,
        hist_range=(0, 3),
        stat="count"
    )


def plot_violin_signed(dataframes, y_ax, out_file='violin.png', what_for='talk'):
    """
    Generate violin plots of the mean signed errors
    of force field energies with respect to QM energies.

    Parameters
    ----------
    msds : 2D numpy array
        Mean signed deviations for each method with reference to first method
        msds[i][j] represents ith mol, jth method's MSE
    y_ax : string
        column of the results df containing the results to plotted in the y axis
    ff_list : list
        list of methods corresponding to energies in msds
    what_for : string
        dictates figure size, text size of axis labels, legend, etc.
        "paper" or "talk"

    """

    # create dataframe from list of lists
    df = pd.DataFrame({method: result.loc[:, y_ax] for method, result in dataframes.items()})

    medians = df.median(axis=0)

    # reshape for grouped violin plots
    df = df.melt(value_vars=list(df.columns))
    df = df.dropna()
    
    # set grid background style
    sns.set(style="whitegrid")

    if what_for == 'paper':
#        f, ax = plt.subplots(figsize=(2, 5))
        f, ax = plt.subplots(figsize=(8, 5))
        large_font = 10
        small_font = 8
        lw = 1
        med_pt = 5
        xrot = 45
        xha = 'right'
    elif what_for == 'talk':
        f, ax = plt.subplots(figsize=(10, 8))
        large_font = 16
        small_font = 14
        lw = 2
        med_pt = 10
        xrot = 0
        xha = 'center'

        # overlapping violins
        #lw=1.0
        #f, ax = plt.subplots(figsize=(4, 8))

    my_cmap = "tab10"
    colors = sns.color_palette(my_cmap)

    # show each distribution with both violins and points
    ax = sns.violinplot(x='variable', y='value', data=df, inner="box",
                        palette=colors, size=5, aspect=3, linewidth=lw)

    # replot the median point for larger marker, zorder to plot points on top
    xlocs = ax.get_xticks()
    for i, x in enumerate(xlocs):
        plt.scatter(x, medians[i], marker='o', color='white', s=med_pt, zorder=100)

    # represent the y-data on log scale
    plt.yscale('symlog')

    # set alpha transparency
    plt.setp(ax.collections, alpha=0.8)

    # hide vertical plot boundaries
    sns.despine(left=True)

    # add labels and adjust font sizes
    ax.set_xlabel("")
    ax.set_ylabel(y_ax, size=large_font)
    plt.xticks(fontsize=small_font, rotation=xrot, ha=xha)

    # settings for overlapping violins
    # plt.xticks([])
    # locs, labels = plt.yticks(fontsize=large_font)
    # plt.yticks(locs, [])
    # plt.xlim(-1, 1)
    # ax.set_ylabel("", size=large_font)

    # save and close figure
    plt.savefig(out_file, dpi=600, bbox_inches='tight')
    #plt.show()
    plt.clf()
    plt.close(plt.gcf())

    # reset plot parameters (white grid)
    sns.reset_orig()


def plot_mol_rmses(mol_name, rmses, xticklabels, eff_nconfs, ref_nconfs, what_for='talk'):
    """
    Generate bar plot of RMSEs of conformer energies for this molecule of
    all methods compared to reference method.

    Number of conformers used to calculate the RMSE is also plotted
    as a solid line. The number of possible conformers available
    by the reference method is plotted as a dashed line.

    Parameters
    ----------
    mol_name : string
        title of the mol being plotted
    rmses : list
        rmses[i] is the RMSE of this mol of reference compared to ith method
    xticklabels : list
        list of methods of the same length as rmses list; should not include
        reference method label
    eff_nconfs : list
        effective number of conformers with non-nan values;
        same format and length as rmses
    ref_nconfs : int
        number of conformers in the reference method
    what_for : string
        dictates figure size, text size of axis labels, legend, etc.
        "paper" or "talk"

    """

    # set figure-related labels
    plttitle = f"RMSEs of relative energies for\nmolecule {mol_name}"
    ylabel = "RMSE (kcal/mol)"
    figname = f"barRMSE_{mol_name}.png"

    # define x locations by integer number of methods
    x_locs = list(range(len(xticklabels)))

    if what_for == 'paper':
        fig = plt.gcf()
        fig.set_size_inches(4, 3)
        large_font = 14
        small_font = 10
    elif what_for == 'talk':
        large_font = 18
        small_font = 16

    # label figure; label xticks before plot for better spacing
    plt.title(plttitle, fontsize=large_font)
    plt.ylabel(ylabel, fontsize=large_font)
    plt.xticks(x_locs, xticklabels, fontsize=small_font, rotation=-30, ha='left')
    plt.yticks(fontsize=small_font)

    # define custom colors here if desired
    #colors = ['tab:blue']*len(xticklabels)
    colors = ['tab:blue']*2 + ['tab:orange']*2 + ['tab:green']*2

    # plot rmses as bars
    plt.bar(x_locs, rmses, color=colors, align='center', label='RMSE')

    # plot number of conformers as lines
    ax2 = plt.twinx()
    ax2.plot(x_locs, eff_nconfs, color='k', alpha=0.5,
             label='actual num confs')
    ax2.axhline(ref_nconfs, color='k', alpha=0.5, ls='--',
                label='reference num confs')

    # format line graph properties, then add plot legend
    ax2.set_ylabel('Number of conformers', fontsize=large_font)
    ax2.tick_params(axis='y', labelsize=small_font)
    ax2.yaxis.set_ticks(np.arange(min(eff_nconfs)-1, ref_nconfs+2, 1))
    plt.legend()

    # save and close figure
    plt.savefig(figname, bbox_inches='tight')
    #plt.show()
    plt.clf()
    plt.close(plt.gcf())


def plot_mol_minima(dataframes, y_label, y_data, method, out_file='minimaE.png', what_for='talk'):
    """
    Generate line plot of conformer energies of all methods (single molecule).

    Parameters
    ----------
    dataframes : Dict[str, pd.DataFrame]
        Dict of Dataframes with method as keys, results as values.
        This function modifies the Dataframes in-place.
    y_label : str
        y-axis label to use on the resulting plot.
    y_data : str
        Column name in each Dataframe corresponding to y-axis data.
    method : str
        Key in `dataframes` used as the reference method for figure generation.
        Can ultimately be any key in the dictionary; does not matter which.
    out_file : Path-like
        Name of output file for rendered figure.
        Extension indicates type of image format used.
    what_for : string
        dictates figure size, text size of axis labels, legend, etc.
        "paper" or "talk"

    """
    if what_for == 'paper':
        fig = plt.figure(figsize=(7.5, 3))
        large_font = 12
        small_font = 10
        xaxis_font = 6
        mark_size = 5
    elif what_for == 'talk':
        fig = plt.figure(figsize=(20, 8))
        large_font = 18
        small_font = 14
        xaxis_font = 10
        mark_size = 9

    for mid in dataframes[method].molecule_index.unique():    
        ref_confs = dataframes[method].loc[dataframes[method].molecule_index==mid]
        ref_nconfs = ref_confs.shape[0]

        # set figure-related labels
        mol_name = list(ref_confs.index)[0]
        plttitle = f"Relative Energies of {mol_name} Minima"
        ylabel = y_label
        figname = f"{out_file[:-4]}_{mol_name}{out_file[-4:]}"
        xlabs = ref_confs['conformer_index']

        # create figure
        ax = fig.gca()
        # label figure; label xticks before plotting for better spacing.
        plt.title(plttitle, fontsize=large_font)
        plt.ylabel(ylabel, fontsize=large_font)
        plt.xlabel("conformer index", fontsize=large_font)
        plt.xticks(list(range(ref_nconfs)), xlabs, fontsize=xaxis_font)
        plt.yticks(fontsize=small_font)

        # define line colors and markers
        colors = mpl.cm.rainbow(np.linspace(0, 1, len(dataframes)))
        markers = [
            "x", "^", "8", "d", "o", "s", "*", "p", "v", "<", "D", "+", ">", "."
        ] * 10
        
        # plot the data
        for i, method in enumerate(dataframes):
            # define x's from integer range with step 1
            xi = list(range(ref_nconfs))
            # get the relative energies of method
            energies = dataframes[method].loc[dataframes[method].molecule_index==mid][y_data]
            # plot values
            plt.plot(xi, energies, color=colors[i], label=method,
                     marker=markers[i], markersize=mark_size, alpha=0.6)

        # add legend, set plot limits, add grid
        plt.legend(bbox_to_anchor=(1, 1), loc=2, prop={'size': small_font})

        plt.grid()

        # save and close figure
        plt.savefig(figname, bbox_inches='tight')
        plt.clf()
        
