# Copyright (C) 2022 National Center for Atmospheric Research and National Oceanic and Atmospheric Administration
# SPDX-License-Identifier: Apache-2.0
#
# Code to create plots for satellite observations
# Copied from surfplots and altered to use xarray syntax instead of pandas


import cartopy.crs as ccrs
import matplotlib as mpl
import matplotlib.pyplot as plt
import monet as monet
import numpy as np
import pandas as pd
import seaborn as sns
import xarray as xr
from numpy import corrcoef

from monet.plots.taylordiagram import TaylorDiagram as td
from monet.util.tools import get_epa_region_bounds as get_epa_bounds
from monet.util.tools import get_giorgi_region_bounds as get_giorgi_bounds

from ..plots import savefig

plt.set_loglevel(level="warning")

sns.set_context("paper")


def time_average(dset, varname=None, period="1D", time_offset=None):
    """Calculates 24-hour averages

    Parameters
    ----------
    dset : dataframe
        Model/obs pair of hourly data
    varname : None | str
        Column label of observation variable to apply the calculation
    period : str
        The period over which to average. Should be noted in Pandas style
        (e.g., '1D', '1h', 'ME', etc.)
    time_offset : None | timedelta
        Whether a time_offset should be applied. Can be useful if accounting
        for time offsets is desirable.

    Returns
    -------
    dataframe
        dataframe with applied calculation
    """
    daily = dset[varname].resample(time=period, offset=time_offset).mean()

    return daily


def map_projection(f):
    """Defines map projection. This needs updating to make it more generic.

    Parameters
    ----------
    f : class
        model class

    Returns
    -------
    cartopy projection
        projection to be used by cartopy in plotting

    """
    import cartopy.crs as ccrs

    if f.model.lower() == "cmaq":
        proj = ccrs.LambertConformal(central_longitude=f.obj.XCENT, central_latitude=f.obj.YCENT)
    elif f.model.lower() == "wrfchem" or f.model.lower() == "rapchem":
        if f.obj.MAP_PROJ == 1:
            proj = ccrs.LambertConformal(
                central_longitude=f.obj.CEN_LON, central_latitude=f.obj.CEN_LAT
            )
        elif f.MAP_PROJ == 6:
            # Plate Carree is the equirectangular or equidistant cylindrical
            proj = ccrs.PlateCarree(central_longitude=f.obj.CEN_LON)
        else:
            raise NotImplementedError(
                "WRFChem projection not supported. Please add to surfplots.py"
            )
    # Need to add the projections you want to use for the other models here.
    elif f.model.lower() == "rrfs":
        proj = ccrs.LambertConformal(
            central_longitude=f.obj.cen_lon, central_latitude=f.obj.cen_lat
        )
    elif f.model.lower() in ["cesm_fv", "cesm_se", "raqms"]:
        proj = ccrs.PlateCarree()
    elif f.model.lower() == "random":
        proj = ccrs.PlateCarree()
    else:  # Let's change this tomorrow to just plot as lambert conformal if nothing provided.
        raise NotImplementedError(
            "Projection not defined for new model. Please add to surfplots.py"
        )
    return proj


def make_timeseries(
    dset,
    varname=None,
    label=None,
    ax=None,
    avg_window="h",
    ylabel=None,
    vmin=None,
    vmax=None,
    domain_type=None,
    domain_name=None,
    plot_dict=None,
    fig_dict=None,
    text_dict=None,
    debug=False,
):
    """Creates timeseries plot.

    Parameters
    ----------
    dset : xr.Dataset
        model/obs pair data to plot
    varname : str
        Variable label of variable to plot
    label : str
        Name of variable to use in plot legend
    ax : ax
        matplotlib ax from previous occurrence so can overlay obs and model
        results on the same plot
    avg_window : rule
        Pandas resampling rule (e.g., 'h', 'D')
    ylabel : str
        Title of y-axis
    vmin : real number
        Min value to use on y-axis
    vmax : real number
        Max value to use on y-axis
    domain_type : str
        Domain type specified in input yaml file
    domain_name : str
        Domain name specified in input yaml file
    plot_dict : dictionary
        Dictionary containing information about plotting for each pair
        (e.g., color, linestyle, markerstyle)
    fig_dict : dictionary
        Dictionary containing information about figure
    text_dict : dictionary
        Dictionary containing information about text
    debug : boolean
        Whether to plot interactively (True) or not (False). Flag for
        submitting jobs to supercomputer turn off interactive mode.

    Returns
    -------
    ax
        matplotlib ax such that driver.py can iterate to overlay multiple models on the
        same plot

    """
    if plot_dict is None:
        plot_dict = {}
    if not debug:
        plt.ioff()
    # First define items for all plots
    # set default text size
    def_text = dict(fontsize=14)
    if text_dict is not None:
        text_kwargs = {**def_text, **text_dict}
    else:
        text_kwargs = def_text
    # set ylabel to column if not specified.
    if ylabel is None:
        ylabel = varname
    if label is not None:
        plot_dict["label"] = label
    if (label is None) and (label not in plot_dict.keys()):
        plot_dict["label"] = varname
    if vmin is not None and vmax is not None:
        plot_dict["ylim"] = [vmin, vmax]
    # scale the fontsize for the x and y labels by the text_kwargs
    plot_dict["fontsize"] = text_kwargs["fontsize"] * 0.8

    # Then, if no plot has been created yet, create a plot and plot the obs.
    if ax is None:
        # First define the colors for the observations.
        obs_dict = dict(color="k", linestyle="-", marker="*", linewidth=1.2, markersize=6.0)
        if plot_dict is not None:
            # Whatever is not defined in the yaml file is filled in with the obs_dict here.
            plot_kwargs = {**obs_dict, **plot_dict}
        else:
            plot_kwargs = obs_dict
        # create the figure
        if fig_dict is not None:
            f, ax = plt.subplots(**fig_dict)
        else:
            f, ax = plt.subplots(figsize=(10, 6))
        # plot the line
        print(plot_kwargs)

        if avg_window is None:
            dset[varname].mean(dim=("y", "x"), skipna=True).plot.line(
                x="time",
                ax=ax,
                color=plot_kwargs["color"],
                linestyle=plot_kwargs["linestyle"],
                marker=plot_kwargs["marker"],
                linewidth=plot_kwargs["linewidth"],
                markersize=plot_kwargs["markersize"],
                label=plot_kwargs["label"],
            )
        else:
            dset[varname].resample(time=avg_window).mean().mean(dim=("y", "x")).plot.line(
                x="time",
                ax=ax,
                color=plot_kwargs["color"],
                linestyle=plot_kwargs["linestyle"],
                marker=plot_kwargs["marker"],
                linewidth=plot_kwargs["linewidth"],
                markersize=plot_kwargs["markersize"],
                label=plot_kwargs["label"],
            )

    # If plot has been created add to the current axes.
    else:
        # this means that an axis handle already exists and use it to plot the model output.
        mod_dict = dict(color=None, linestyle="-", marker="*", linewidth=1.2, markersize=6.0)
        if plot_dict is not None:
            # Whatever is not defined in the yaml file is filled in with the mod_dict here.
            plot_kwargs = {**mod_dict, **plot_dict}
        else:
            plot_kwargs = obs_dict
        if avg_window is None:
            dset[varname].mean(dim=("y", "x")).plot.line(
                x="time",
                ax=ax,
                color=plot_kwargs["color"],
                linestyle=plot_kwargs["linestyle"],
                marker=plot_kwargs["marker"],
                linewidth=plot_kwargs["linewidth"],
                markersize=plot_kwargs["markersize"],
                label=plot_kwargs["label"],
            )
        else:
            dset[varname].resample(time=avg_window).mean().mean(dim=("y", "x")).plot.line(
                x="time",
                ax=ax,
                color=plot_kwargs["color"],
                linestyle=plot_kwargs["linestyle"],
                marker=plot_kwargs["marker"],
                linewidth=plot_kwargs["linewidth"],
                markersize=plot_kwargs["markersize"],
                label=plot_kwargs["label"],
            )

    # Set parameters for all plots
    ax.set_ylabel(ylabel, fontweight="bold", **text_kwargs)
    ax.set_xlabel("time", fontweight="bold", **text_kwargs)
    ax.legend(frameon=False, fontsize=text_kwargs["fontsize"] * 0.8)
    ax.tick_params(axis="both", length=10.0, direction="inout")
    ax.tick_params(axis="both", which="minor", length=5.0, direction="out")
    ax.legend(
        frameon=False,
        fontsize=text_kwargs["fontsize"] * 0.8,
        bbox_to_anchor=(1.0, 0.9),
        loc="center left",
    )
    if domain_type is not None and domain_name is not None:
        if domain_type == "epa_region":
            ax.set_title("EPA Region " + domain_name, fontweight="bold", **text_kwargs)
        else:
            ax.set_title(domain_name, fontweight="bold", **text_kwargs)
    return ax


def make_taylor(
    dset,
    varname_o=None,
    label_o="Obs",
    varname_m=None,
    label_m="Model",
    mean_criteria=None,
    dia=None,
    ylabel=None,
    ty_scale=1.5,
    domain_type=None,
    domain_name=None,
    plot_dict=None,
    fig_dict=None,
    text_dict=None,
    debug=False,
    normalize=False,
    scale_factor=1,
):
    """Creates taylor plot. Note sometimes model values are off the scale
    on this plot. This will be fixed soon.

    Parameters
    ----------
    dset : xr.Dataset
        model/obs pair data to plot
    column_o : str
        Column label of observational variable to plot
    label_o : str
        Name of observational variable to use in plot legend
    column_m : str
        Column label of model variable to plot
    label_m : str
        Name of model variable to use in plot legend
    mean_criteria : str
        'None', 'space', 'time'. If None, values and correlations are compared
        over all dimensions (x, y and time). If 'space', the spatial mean over the
        comain is calculated before doing the comparison. If 'time', the temporal
        mean is calculated before doing the comparison.
    dia : dia
        matplotlib ax from previous occurrence so can overlay obs and model
        results on the same plot
    ylabel : str
        Title of x-axis
    ty_scale : real
        Scale to apply to taylor plot to control the plotting range
    domain_type : str
        Domain type specified in input yaml file
    domain_name : str
        Domain name specified in input yaml file
    plot_dict : dictionary
        Dictionary containing information about plotting for each pair
        (e.g., color, linestyle, markerstyle)
    fig_dict : dictionary
        Dictionary containing information about figure
    text_dict : dictionary
        Dictionary containing information about text
    debug : boolean
        Whether to plot interactively (True) or not (False). Flag for
        submitting jobs to supercomputer turn off interactive mode.

    Returns
    -------
    class
        Taylor diagram class defined in MONET

    """
    # import pdb; pdb.set_trace()

    if mean_criteria == "space":
        dset_forplot = dset.mean(dim=("x", "y"))
    elif mean_criteria == "time":
        dset_forplot = dset.mean(dim="time")
    else:
        dset_forplot = dset
    # import pdb; pdb.set_trace()

    nan_ind = (~np.isnan(dset_forplot[varname_o].values)) & (
        ~np.isnan(dset_forplot[varname_m].values)
    )
    # First define items for all plots
    if not debug:
        plt.ioff()

    # set default text size
    def_text = dict(fontsize=14.0)
    if text_dict is not None:
        text_kwargs = {**def_text, **text_dict}
    else:
        text_kwargs = def_text
    # set ylabel to column if not specified.
    if ylabel is None:
        ylabel = varname_o
    # Then, if no plot has been created yet, create a plot and plot the first pair.
    refstd = dset_forplot[varname_o].std().values

    if dia is None:
        # create the figure
        if fig_dict is not None:
            f = plt.figure(**fig_dict)
        else:
            f = plt.figure(figsize=(12, 10))
        sns.set_style("ticks")
        # plot the line
        cc = corrcoef(
            dset_forplot[varname_o].values[nan_ind].flatten(),
            dset_forplot[varname_m].values[nan_ind].flatten(),
        )[0, 1]
        if normalize:
            print(f"Base standard deviation: {refstd: 1.3g}")
            scale_factor = refstd
            dia = td(1, scale=ty_scale, fig=f, rect=111, label=label_o)
            dia.add_sample(
                dset_forplot[varname_m].std().values / scale_factor,
                cc,
                zorder=9,
                label=label_m,
                **plot_dict,
            )
        elif scale_factor != 1:
            dset_forplot[varname_m].attrs[
                "units"
            ] = f"{scale_factor} {dset_forplot[varname_m].attrs['units']}"

        else:
            dia = td(refstd, scale=ty_scale, fig=f, rect=111, label=label_o)
            dia.add_sample(
                dset_forplot[varname_m].std().values,
                cc,
                zorder=9,
                label=label_m,
                **plot_dict,
            )
        plt.grid(linewidth=1, alpha=0.5)

    # If plot has been created add to the current axes.
    else:
        # this means that an axis handle already exists and use it to plot another model
        cc = corrcoef(
            dset_forplot[varname_o].values[nan_ind].flatten(),
            dset_forplot[varname_m].values[nan_ind].flatten(),
        )[0, 1]
        if normalize:
            scale_factor = refstd
            dia.add_sample(
                dset_forplot[varname_m].std().values / scale_factor,
                cc,
                zorder=9,
                label=label_m,
                **plot_dict,
            )

        else:
            dia.add_sample(
                dset_forplot[varname_m].std().values,
                cc,
                zorder=9,
                label=label_m,
                **plot_dict,
            )
    # Set parameters for all plots
    contours = dia.add_contours(colors="0.5")
    # control the clabel format for very high values (e.g., NO2 columns), M.Li
    # plt.clabel(contours, inline=1, fontsize=text_kwargs['fontsize']*0.8)
    plt.clabel(contours, inline=1, fontsize=text_kwargs["fontsize"] * 0.8, fmt="(%1.5g)")

    plt.grid(alpha=0.5)
    plt.legend(
        frameon=False,
        fontsize=text_kwargs["fontsize"] * 0.8,
        bbox_to_anchor=(0.75, 0.93),
        loc="center left",
    )
    if domain_type is not None and domain_name is not None:
        if domain_type == "epa_region":
            plt.title("EPA Region " + domain_name, fontweight="bold", **text_kwargs)
        else:
            plt.title(domain_name, fontweight="bold", **text_kwargs)

    ax = plt.gca()
    ax.axis["left"].label.set_text("Standard Deviation: " + ylabel)
    ax.axis["top"].label.set_text("Correlation")
    ax.axis["left"].label.set_fontsize(text_kwargs["fontsize"])
    ax.axis["top"].label.set_fontsize(text_kwargs["fontsize"])
    ax.axis["left"].label.set_fontweight("bold")
    ax.axis["top"].label.set_fontweight("bold")
    ax.axis["top"].major_ticklabels.set_fontsize(text_kwargs["fontsize"] * 0.8)
    ax.axis["left"].major_ticklabels.set_fontsize(text_kwargs["fontsize"] * 0.8)
    ax.axis["right"].major_ticklabels.set_fontsize(text_kwargs["fontsize"] * 0.8)
    return dia


def calculate_boxplot(
    dset,
    varname=None,
    label=None,
    plot_dict=None,
    comb_bx=None,
    label_bx=None,
):
    """Combines data into acceptable format for box-plot

    Parameters
    ----------
    dset : xr.Dataset
         model/obs pair data to plot
    varname : str
        Dataset label of variable to plot
    label : str
        Name of variable to use in plot legend
    comb_bx: dataframe
        dataframe containing information to create box-plot from previous
        occurrence so can overlay multiple model results on plot
    label_bx: list
        list of string labels to use in box-plot from previous occurrence so
        can overlay multiple model results on plot
    Returns
    -------
    dataframe, list
        dataframe containing information to create box-plot
        list of string labels to use in box-plot

    """
    if comb_bx is None and label_bx is None:
        comb_bx = pd.DataFrame()
        label_bx = []
        # First define the colors for the observations.
        obs_dict = dict(color="gray", linestyle="-", marker="x", linewidth=1.2, markersize=6.0)
        if plot_dict is not None:
            # Whatever is not defined in the yaml file is filled in with the obs_dict here.
            plot_kwargs = {**obs_dict, **plot_dict}
        else:
            plot_kwargs = obs_dict
    else:
        plot_kwargs = plot_dict
    # For all, a column to the dataframe and append the label info to the list.
    plot_kwargs["varname"] = varname
    plot_kwargs["label"] = label
    comb_bx[label] = dset[varname].mean(dim=("x", "y")).to_dataframe()
    label_bx.append(plot_kwargs)

    return comb_bx, label_bx


def make_boxplot(
    comb_bx,
    label_bx,
    ylabel=None,
    vmin=None,
    vmax=None,
    outname="plot",
    domain_type=None,
    domain_name=None,
    plot_dict=None,
    fig_dict=None,
    text_dict=None,
    debug=False,
):
    """Creates box-plot.

    Parameters
    ----------
    comb_bx: dataframe
        dataset containing information to create box-plot from
        calculate_boxplot
    label_bx: list
        list of string labels to use in box-plot from calculate_boxplot
    ylabel : str
        Title of y-axis
    vmin : real number
        Min value to use on y-axis
    vmax : real number
        Max value to use on y-axis
    outname : str
        file location and name of plot (do not include .png)
    domain_type : str
        Domain type specified in input yaml file
    domain_name : str
        Domain name specified in input yaml file
    plot_dict : dictionary
        Dictionary containing information about plotting for each pair
        (e.g., color, linestyle, markerstyle)
    fig_dict : dictionary
        Dictionary containing information about figure
    text_dict : dictionary
        Dictionary containing information about text
    debug : boolean
        Whether to plot interactively (True) or not (False). Flag for
        submitting jobs to supercomputer turn off interactive mode.

    Returns
    -------
    plot
        box plot

    """
    if not debug:
        plt.ioff()
    # First define items for all plots
    # set default text size
    def_text = dict(fontsize=14)
    if text_dict is not None:
        text_kwargs = {**def_text, **text_dict}
    else:
        text_kwargs = def_text
    # set ylabel to column if not specified.
    if ylabel is None:
        ylabel = label_bx[0]["label"]

    # Fix the order and palate colors
    order_box = []
    pal = {}
    for i in range(len(label_bx)):
        order_box.append(label_bx[i]["label"])
        pal[label_bx[i]["label"]] = label_bx[i]["color"]

    # Make plot
    if fig_dict is not None:
        f, ax = plt.subplots(**fig_dict)
    else:
        f, ax = plt.subplots(figsize=(8, 8))
    # Define characteristics of boxplot.
    boxprops = {"edgecolor": "k", "linewidth": 1.5}
    lineprops = {"color": "k", "linewidth": 1.5}
    boxplot_kwargs = {
        "boxprops": boxprops,
        "medianprops": lineprops,
        "whiskerprops": lineprops,
        "capprops": lineprops,
        "fliersize": 2.0,
        "flierprops": dict(
            marker="*", markerfacecolor="blue", markeredgecolor="none", markersize=6.0
        ),
        "width": 0.75,
        "palette": pal,
        "order": order_box,
        "showmeans": True,
        "meanprops": {
            "marker": ".",
            "markerfacecolor": "black",
            "markeredgecolor": "black",
            "markersize": 20.0,
        },
    }
    sns.set_style("whitegrid")
    sns.set_style("ticks")
    sns.boxplot(
        ax=ax,
        x="variable",
        y="value",
        hue="variable",
        data=pd.melt(comb_bx),
        **boxplot_kwargs,
    )
    ax.set_xlabel("")
    ax.set_ylabel(ylabel, fontweight="bold", **text_kwargs)
    ax.tick_params(labelsize=text_kwargs["fontsize"] * 0.8)
    ax.yaxis.get_offset_text().set_fontsize(text_kwargs["fontsize"] * 0.8)
    if domain_type is not None and domain_name is not None:
        if domain_type == "epa_region":
            ax.set_title("EPA Region " + domain_name, fontweight="bold", **text_kwargs)
        else:
            ax.set_title(domain_name, fontweight="bold", **text_kwargs)
    if vmin is not None and vmax is not None:
        ax.set_ylim(ymin=vmin, ymax=vmax)

    plt.tight_layout()
    savefig(
        outname + ".png",
        loc=4,
        logo_height=100,
        decorate=True,
        bbox_inches="tight",
        dpi=200,
    )


def make_spatial_dist(
    dset,
    varname=None,
    label=None,
    ylabel=None,
    vmin=None,
    vmax=None,
    nlevels=None,
    proj=None,
    outname="plot",
    domain_type=None,
    domain_name=None,
    fig_dict=None,
    text_dict=None,
    debug=False,
):
    """Creates difference plot for satellite and model data.
    For data in swath format, overplots all differences
    For data on regular grid, mean difference.

    Parameters
    ----------
    dset : xr.Dataset
        Dataset containing the paired data

    """
    if not debug:
        plt.ioff()

    def_map = dict(states=True, figsize=[15, 8])
    if fig_dict is not None:
        map_kwargs = {**def_map, **fig_dict}
    else:
        map_kwargs = def_map

    # set default text size
    def_text = dict(fontsize=20)
    if text_dict is not None:
        text_kwargs = {**def_text, **text_dict}
    else:
        text_kwargs = def_text

    # set ylabel to column if not specified.
    if ylabel is None:
        ylabel = varname

    # Take the difference for the model output - the sat output

    var2plot = dset[varname]  # Take mean over time,

    if len(var2plot.dims) == 3:
        var2plot = var2plot.mean("time")

    # Determine the domain
    if domain_type == "all" and domain_name == "CONUS":
        latmin = 25.0
        lonmin = -130.0
        latmax = 50.0
        lonmax = -60.0
        title_add = domain_name + ": "
    elif domain_type == "epa_region" and domain_name is not None:
        latmin, lonmin, latmax, lonmax, acro = get_epa_bounds(index=None, acronym=domain_name)
        title_add = "EPA Region " + domain_name + ": "
    elif domain_type == "giorgi_region" and domain_name is not None:
        latmin, lonmin, latmax, lonmax, acro = get_giorgi_bounds(index=None, acronym=domain_name)
        title_add = "Giorgi Region " + domain_name + ": "
    else:
        latmin = -90
        lonmin = -180
        latmax = 90
        lonmax = 180
        title_add = domain_name + ": "

    # Map the model output first.
    cbar_kwargs = dict(aspect=15, shrink=0.8)

    # Add options that this could be included in the fig_kwargs in yaml file too.
    if "extent" not in map_kwargs:
        map_kwargs["extent"] = [lonmin, lonmax, latmin, latmax]
    if "crs" not in map_kwargs:
        map_kwargs["crs"] = proj

    # First determine colorbar
    if vmin is None and vmax is None:
        # vmin = vmodel_mean.quantile(0.01)
        vmax = np.max(
            (
                np.abs(var2plot.quantile(0.99)),
                np.abs(var2plot.quantile(0.01)),
            )
        )
        vmin = -vmax

    if nlevels is None:
        nlevels = 21
    print(vmin, vmax)
    clevel = np.linspace(vmin, vmax, nlevels)
    cmap = mpl.cm.get_cmap("bwr", nlevels - 1)
    norm = mpl.colors.BoundaryNorm(clevel, ncolors=cmap.N, clip=False)

    # I add extend='both' here because the colorbar is setup to plot the values outside the range
    ax = monet.plots.mapgen.draw_map(
        crs=map_kwargs["crs"], extent=map_kwargs["extent"], states=True, counties=True
    )
    # draw scatter plot of model and satellite differences
    c = ax.axes.scatter(dset.longitude, dset.latitude, c=var2plot, cmap=cmap, s=2, norm=norm)
    plt.gcf().canvas.draw()
    plt.tight_layout(pad=0)
    plt.title(title_add + label, fontweight="bold", **text_kwargs)
    ax.axes.set_extent(map_kwargs["extent"], crs=ccrs.PlateCarree())

    # Uncomment these lines if you update above just to verify colorbars are identical.
    # Also specify plot above scatter = ax.axes.scatter etc.
    # cbar = ax.figure.get_axes()[1]
    plt.colorbar(c, ax=ax, extend="both")

    # Update colorbar
    f = plt.gcf()

    model_ax = f.get_axes()[0]
    cax = f.get_axes()[1]

    # get the position of the plot axis and use this to rescale nicely the color bar to the height of the plot.
    position_m = model_ax.get_position()
    position_c = cax.get_position()
    cax.set_position(
        [
            position_c.x0,
            position_m.y0,
            position_c.x1 - position_c.x0,
            (position_m.y1 - position_m.y0) * 1.1,
        ]
    )
    cax.set_ylabel(ylabel, fontweight="bold", **text_kwargs)
    cax.tick_params(
        labelsize=text_kwargs["fontsize"] * 0.8,
        length=10.0,
        width=2.0,
        grid_linewidth=2.0,
    )

    cax.yaxis.get_offset_text().set_fontsize(text_kwargs["fontsize"] * 0.8)

    # plt.tight_layout(pad=0)
    savefig(
        outname + ".png",
        loc=4,
        logo_height=100,
        decorate=True,
        bbox_inches="tight",
        dpi=150,
    )
    return ax


def make_spatial_bias_gridded(
    dset,
    varname_o=None,
    label_o=None,
    varname_m=None,
    label_m=None,
    ylabel=None,
    vmin=None,
    vmax=None,
    nlevels=None,
    proj=None,
    outname="plot",
    domain_type=None,
    domain_name=None,
    fig_dict=None,
    text_dict=None,
    debug=False,
):
    """Creates difference plot for satellite and model data.
    For data in swath format, overplots all differences
    For data on regular grid, mean difference.

    Parameters
    ----------
    dset : xr.Dataset
        Dataset containing the paired data

    """
    if not debug:
        plt.ioff()

    def_map = dict(states=True, figsize=[15, 8])
    if fig_dict is not None:
        map_kwargs = {**def_map, **fig_dict}
    else:
        map_kwargs = def_map

    # set default text size
    def_text = dict(fontsize=20)
    if text_dict is not None:
        text_kwargs = {**def_text, **text_dict}
    else:
        text_kwargs = def_text

    # set ylabel to column if not specified.
    if ylabel is None:
        ylabel = varname_o

    # Take the difference for the model output - the sat output

    diff_mod_min_obs = (dset[varname_m] - dset[varname_o]).squeeze()
    # Take mean over time,
    if len(diff_mod_min_obs.dims) == 3:
        diff_mod_min_obs = diff_mod_min_obs.mean("time")

    # Determine the domain
    if domain_type == "all" and domain_name == "CONUS":
        latmin = 25.0
        lonmin = -130.0
        latmax = 50.0
        lonmax = -60.0
        title_add = domain_name + ": "
    elif domain_type == "epa_region" and domain_name is not None:
        latmin, lonmin, latmax, lonmax, acro = get_epa_bounds(index=None, acronym=domain_name)
        title_add = "EPA Region " + domain_name + ": "
    elif domain_type == "giorgi_region" and domain_name is not None:
        latmin, lonmin, latmax, lonmax, acro = get_giorgi_bounds(index=None, acronym=domain_name)
        title_add = "Giorgi Region " + domain_name + ": "
    else:
        latmin = -90
        lonmin = -180
        latmax = 90
        lonmax = 180
        title_add = domain_name + ": "

    # Map the model output first.
    cbar_kwargs = dict(aspect=15, shrink=0.8)

    # Add options that this could be included in the fig_kwargs in yaml file too.
    if "extent" not in map_kwargs:
        map_kwargs["extent"] = [lonmin, lonmax, latmin, latmax]
    if "crs" not in map_kwargs:
        map_kwargs["crs"] = proj

    # First determine colorbar
    if vmin is None and vmax is None:
        # vmin = vmodel_mean.quantile(0.01)
        vmax = np.max(
            (
                np.abs(diff_mod_min_obs.quantile(0.99)),
                np.abs(diff_mod_min_obs.quantile(0.01)),
            )
        )
        vmin = -vmax

    if nlevels is None:
        nlevels = 21
    print(vmin, vmax)
    clevel = np.linspace(vmin, vmax, nlevels)
    cmap = mpl.cm.get_cmap("bwr", nlevels - 1)
    norm = mpl.colors.BoundaryNorm(clevel, ncolors=cmap.N, clip=False)

    # I add extend='both' here because the colorbar is setup to plot the values outside the range
    ax = monet.plots.mapgen.draw_map(
        crs=map_kwargs["crs"], extent=map_kwargs["extent"], states=True, counties=True
    )
    # draw scatter plot of model and satellite differences
    c = ax.axes.scatter(
        dset.longitude, dset.latitude, c=diff_mod_min_obs, cmap=cmap, s=2, norm=norm
    )
    plt.gcf().canvas.draw()
    plt.tight_layout(pad=0)
    plt.title(title_add + label_m + " - " + label_o, fontweight="bold", **text_kwargs)
    ax.axes.set_extent(map_kwargs["extent"], crs=ccrs.PlateCarree())

    # Uncomment these lines if you update above just to verify colorbars are identical.
    # Also specify plot above scatter = ax.axes.scatter etc.
    # cbar = ax.figure.get_axes()[1]
    plt.colorbar(c, ax=ax, extend="both")

    # Update colorbar
    f = plt.gcf()

    model_ax = f.get_axes()[0]
    cax = f.get_axes()[1]

    # get the position of the plot axis and use this to rescale nicely the color bar to the height of the plot.
    position_m = model_ax.get_position()
    position_c = cax.get_position()
    cax.set_position(
        [
            position_c.x0,
            position_m.y0,
            position_c.x1 - position_c.x0,
            (position_m.y1 - position_m.y0) * 1.1,
        ]
    )
    cax.set_ylabel("$\Delta$" + ylabel, fontweight="bold", **text_kwargs)
    cax.tick_params(
        labelsize=text_kwargs["fontsize"] * 0.8,
        length=10.0,
        width=2.0,
        grid_linewidth=2.0,
    )

    cax.yaxis.get_offset_text().set_fontsize(text_kwargs["fontsize"] * 0.8)

    # plt.tight_layout(pad=0)
    savefig(
        outname + ".png",
        loc=4,
        logo_height=100,
        decorate=True,
        bbox_inches="tight",
        dpi=150,
    )
    return ax


def make_multi_boxplot(
    comb_bx,
    label_bx,
    region_bx,
    region_list=None,
    model_name_list=None,
    ylabel=None,
    vmin=None,
    vmax=None,
    outname="plot",
    domain_type=None,
    domain_name=None,
    plot_dict=None,
    fig_dict=None,
    text_dict=None,
    debug=False,
):
    """Creates box-plot.

    Parameters
    ----------
    comb_bx : dataframe
        dataframe containing information to create box-plot from
        calculate_boxplot
    label_bx : list
        list of string labels to use in box-plot from calculate_boxplot
    region_bx : dataframe
        dataframe containing information of boxes to help create multi-box-plot
        from calculate_boxplot
    model_name_list : list of str
        list of models and observation sources used for x-labels in plot
    ylabel : str
        Title of y-axis
    vmin : real number
        Min value to use on y-axis
    vmax : real number
        Max value to use on y-axis
    outname : str
        file location and name of plot (do not include .png)
    domain_type : str
        Domain type specified in input yaml file
    domain_name : str
        Domain name specified in input yaml file
    plot_dict : dictionary
        Dictionary containing information about plotting for each pair
        (e.g., color, linestyle, markerstyle)
    fig_dict : dictionary
        Dictionary containing information about figure
    text_dict : dictionary
        Dictionary containing information about text
    debug : boolean
        Whether to plot interactively (True) or not (False). Flag for
        submitting jobs to supercomputer turn off interactive mode.

    Returns
    -------
    plot
        multi-box plot

    """
    if not debug:
        plt.ioff()
    # First define items for all plots
    # set default text size
    def_text = dict(fontsize=14)
    if text_dict is not None:
        text_kwargs = {**def_text, **text_dict}
    else:
        text_kwargs = def_text
    # set ylabel to column if not specified.
    if ylabel is None:
        ylabel = label_bx[0]["column"]

    # Fix the order and palate colors
    order_box = []
    pal = {}
    for i in range(len(label_bx)):
        order_box.append(label_bx[i]["label"])
        pal[label_bx[i]["label"]] = label_bx[i]["color"]

    # Make plot
    if fig_dict is not None:
        f, ax = plt.subplots(**fig_dict)
    else:
        f, ax = plt.subplots(figsize=(8, 8))
    # Define characteristics of boxplot.
    sns.set_style("whitegrid")
    sns.set_style("ticks")
    len_combx = len(comb_bx.columns)

    data_obs = comb_bx[comb_bx.columns[0]].to_frame().rename({comb_bx.columns[0]: "Value"}, axis=1)
    data_obs["model"] = model_name_list[0]
    data_obs["Regions"] = region_bx["set_regions"].values
    to_concat = []
    to_concat.append(data_obs[["Value", "model", "Regions"]])

    for i in range(1, len_combx):
        data_model = (
            comb_bx[comb_bx.columns[i]].to_frame().rename({comb_bx.columns[i]: "Value"}, axis=1)
        )
        data_model["model"] = model_name_list[i]
        data_model["Regions"] = region_bx["set_regions"].values
        to_concat.append(data_model[["Value", "model", "Regions"]])

    tdf = pd.concat(to_concat)
    acro = region_list
    sns.boxplot(
        x="Regions",
        y="Value",
        hue="model",
        data=tdf.loc[tdf.Regions.isin(acro)],
        order=acro,
        showfliers=False,
    )
    ax.set_xlabel("")
    ax.set_ylabel(ylabel, fontweight="bold", **text_kwargs)
    ax.tick_params(labelsize=text_kwargs["fontsize"] * 0.8)
    if domain_type is not None and domain_name is not None:
        if domain_type == "epa_region":
            ax.set_title("EPA Region " + domain_name, fontweight="bold", **text_kwargs)
        elif domain_type == "giorgi_region":
            ax.set_title("Giorgi Region" + domain_name, fontweight="bold", **text_kwargs)
        else:
            ax.set_title(domain_name, fontweight="bold", **text_kwargs)
    if vmin is not None and vmax is not None:
        ax.set_ylim(ymin=vmin, ymax=vmax)

    plt.tight_layout()
    savefig(outname + ".png", loc=4, logo_height=100)
