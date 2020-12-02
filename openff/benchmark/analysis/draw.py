def draw_scatter(
    x_data, y_data, method_labels, x_label, y_label, out_file, what_for="talk"
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
    print(f"\nNumber of data points in full scatter plot: {len(flatten(x_data))}")
    markers = ["o", "^", "d", "x", "s", "p", "P", "3", ">"]

    num_methods = len(x_data)
    plist = []
    for i in range(num_methods):
        p = plt.scatter(
            x_data[i],
            y_data[i],
            marker=markers[i],
            label=method_labels[i + 1],
            alpha=0.6,
        )
        plist.append(p)

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
        for p in plist:
            p.set_sizes([8.0])

    elif what_for == "talk":
        plt.xlabel(x_label, fontsize=14)
        plt.ylabel(y_label, fontsize=14)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.legend(loc=(1.04, 0.5), fontsize=12)
        # make the marker size smaller
        for p in plist:
            p.set_sizes([4.0])

    # set log scaling but use symmetric log for negative values
    #    plt.yscale('symlog')

    plt.savefig(out_file, bbox_inches="tight")
    plt.clf()
    # plt.show()





def draw_corr(
    x_data, y_data, method_labels, x_label, y_label, out_file, what_for="talk"
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
    print(f"\nNumber of data points in full scatter plot: {len(flatten(x_data))}")
    markers = ["o", "^", "d", "x", "s", "p", "P", "3", ">"]

    num_methods = len(x_data)
    plist = []
    for i in range(num_methods):
        p = plt.scatter(
            x_data[i],
            y_data[i],
            marker=markers[i],
            label=method_labels[i + 1],
            alpha=0.6,
        )
        plist.append(p)

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
        for p in plist:
            p.set_sizes([8.0])

    elif what_for == "talk":
        plt.xlabel(x_label, fontsize=14)
        plt.ylabel(y_label, fontsize=14)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.legend(loc=(1.04, 0.5), fontsize=12)
        # make the marker size smaller
        for p in plist:
            p.set_sizes([4.0])

    plt.savefig(out_file, bbox_inches="tight")
    plt.clf()
    # plt.show()


def draw_ridgeplot(
    mydata,
    method_labels,
    x_label,
    out_file,
    what_for="paper",
    bw="scott",
    same_subplot=False,
    sym_log=False,
    hist_range=(-15, 15),
):

    """
    Draw ridge plot of data (to which kernel density estimate is applied)
    segregated by each method (representing a different color/level).

    Modified from the following code:
    https://seaborn.pydata.org/examples/kde_ridgeplot.html

    Parameters
    ----------
    mydata : list of lists
        mydata[i][j] represents ith method and jth molecular structure
    method_labels : list
        list of all the method names including reference method first
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

    num_methods = len(mydata)

    # Initialize the FacetGrid object
    my_cmap = "tab10"
    sns.palplot(sns.color_palette(my_cmap))
    colors = sns.color_palette(my_cmap)

    # convert data to dataframes for ridge plot
    temp = []
    for i in range(num_methods):
        df = pd.DataFrame(mydata[i], columns=[x_label])
        df["method"] = method_labels[i + 1]
        temp.append(df)

    # list of dataframes concatenated to single dataframe
    df = pd.concat(temp, ignore_index=True)
    #    print(method_labels)
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
                sns.distplot,
                x_label,
                hist=True,
                kde=False,
                bins=15,
                norm_hist=True,
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
            "histtype": "step",
            "alpha": 1.0,
            "linewidth": ridgedict["lw"],
            "range": hist_range,
            "align": "mid",
        }
        g.map(
            sns.distplot,
            x_label,
            hist=True,
            kde=False,
            bins=15,
            norm_hist=True,
            hist_kws=histoptions,
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
        n_ffs = len(method_labels) - 1
        for i in range(n_ffs):
            patches.append(
                mpl.patches.Patch(
                    color=cmap(i/10),
                    label=method_labels[i + 1],
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
    plt.ylabel("Density", fontsize=ridgedict["xfontsize"])
    plt.xticks(fontsize=ridgedict["xfontsize"])

    # save with transparency for overlapping plots
    plt.savefig(out_file, transparent=True, bbox_inches="tight")
    plt.clf()
    # plt.show()


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
