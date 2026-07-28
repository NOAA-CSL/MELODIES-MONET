# SPDX-License-Identifier: Apache-2.0
#
#Code to create plots for satellite observations
# Copied from surfplots and altered to use xarray syntax instead of pandas

import monet as monet
import seaborn as sns
from monet.util.tools import calc_8hr_rolling_max, calc_24hr_ave
import pandas as pd
import numpy as np
import cartopy.crs as ccrs
import matplotlib as mpl
import matplotlib.pyplot as plt
from numpy import corrcoef
sns.set_context('paper')
from monet.plots.taylordiagram import TaylorDiagram as td
from matplotlib.colors import ListedColormap
from monet.util.tools import get_epa_region_bounds as get_epa_bounds 
import math
try:
    import uxarray as ux
except ImportError:  # optional dep
    ux = None
from melodies_monet.plots import savefig

def make_24hr_regulatory(df, col=None):
    """Calculates 24-hour averages
    
    Parameters
    ----------
    df : dataframe
        Model/obs pair of hourly data
    col : str
        Column label of observation variable to apply the calculation 
    Returns
    -------
    dataframe
        dataframe with applied calculation
        
    """
    return calc_24hr_ave(df, col)


def make_8hr_regulatory(df, col=None):
    """Calculates 8-hour rolling average daily
    
    Parameters
    ----------
    df : dataframe
        Model/obs pair of hourly data
    col : str
        Column label of observation variable to apply the calculation 
    Returns
    -------
    dataframe
        dataframe with applied calculation
        
    """
    return calc_8hr_rolling_max(df, col, window=8)

def calc_default_colors(p_index):
    """List of default colors, lines, and markers to use if user does not 
    specify them in the input yaml file.
    
    Parameters
    ----------
    p_index : integer
        Number of pairs in analysis class
    
    Returns
    -------
    list
        List of dictionaries containing default colors, lines, and 
        markers to use for plotting for the number of pairs in analysis class
        
    """
    x = [dict(color='b', linestyle='--',marker='x'),
         dict(color='g', linestyle='-.',marker='o'),
         dict(color='r', linestyle=':',marker='v'),
         dict(color='c', linestyle='--',marker='^'),
         dict(color='m', linestyle='-.',marker='s')]
    #Repeat these 5 instances over and over if more than 5 lines.
    return x[p_index % 5]

def new_color_map():
    """Creates new color map for difference plots
    
    Returns
    -------
    colormap
        Orange and blue color map
        
    """
    top = mpl.cm.get_cmap('Blues_r', 128)
    bottom = mpl.cm.get_cmap('Oranges', 128)
    newcolors = np.vstack((top(np.linspace(0, 1, 128)),
                           bottom(np.linspace(0, 1, 128))))
    return ListedColormap(newcolors, name='OrangeBlue')

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
    if f.model.lower() == 'cmaq':
        proj = ccrs.LambertConformal(
            central_longitude=f.obj.XCENT, central_latitude=f.obj.YCENT)
    elif f.model.lower() == 'wrfchem' or f.model.lower() == 'rapchem':
        if f.obj.MAP_PROJ == 1:
            proj = ccrs.LambertConformal(
                central_longitude=f.obj.CEN_LON, central_latitude=f.obj.CEN_LAT)
        elif f.MAP_PROJ == 6:
            #Plate Carree is the equirectangular or equidistant cylindrical
            proj = ccrs.PlateCarree(
                central_longitude=f.obj.CEN_LON)
        else:
            raise NotImplementedError('WRFChem projection not supported. Please add to surfplots.py')         
    #Need to add the projections you want to use for the other models here.        
    elif f.model.lower() in ('rrfs', 'ufs'):
        proj = ccrs.LambertConformal(
            central_longitude=f.obj.cen_lon, central_latitude=f.obj.cen_lat)
    elif f.model.lower() in ['cesm_fv','cesm_se','raqms']:
        proj = ccrs.PlateCarree()
    elif f.model.lower() == 'random':
        proj = ccrs.PlateCarree()
    else: #Let's change this tomorrow to just plot as lambert conformal if nothing provided.
        raise NotImplementedError('Projection not defined for new model. Please add to surfplots.py')
    return proj
    
def make_timeseries(df, df_reg=None,column=None, label=None, ax=None, avg_window=None, ylabel=None,
                    vmin = None, vmax = None,
                    domain_type=None, domain_name=None,
                    plot_dict=None, fig_dict=None, text_dict=None,debug=False):
    """Creates timeseries plot. 
    
    Parameters
    ----------
    df : dataframe
        model/obs pair data to plot
    column : str
        Column label of variable to plot
    df_reg: not currently enabled. empty argument for symmetry with surfplots
        model/obs paired regulatory data to plot
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
    if debug is False:
        plt.ioff()
    #First define items for all plots
    #set default text size
    def_text = dict(fontsize=14)
    if text_dict is not None:
        text_kwargs = {**def_text, **text_dict}
    else:
        text_kwargs = def_text
    # set ylabel to column if not specified.
    if ylabel is None:
        ylabel = column
    if label is not None:
        plot_dict['label'] = label
    if vmin is not None and vmax is not None:
        plot_dict['ylim'] = [vmin,vmax]
    #scale the fontsize for the x and y labels by the text_kwargs
    plot_dict['fontsize'] = text_kwargs['fontsize']*0.8
    
    #Then, if no plot has been created yet, create a plot and plot the obs.
    if ax is None: 
        #First define the colors for the observations.
        obs_dict = dict(color='k', linestyle='-',marker='*', linewidth=1.2, markersize=6.)
        if plot_dict is not None:
            #Whatever is not defined in the yaml file is filled in with the obs_dict here.
            plot_kwargs = {**obs_dict, **plot_dict}
        else:
            plot_kwargs = obs_dict
        # create the figure
        if fig_dict is not None:
            f,ax = plt.subplots(**fig_dict)    
        else: 
            f,ax = plt.subplots(figsize=(10,6))
        # plot the line
        print(plot_kwargs)
        # {'color': 'k', 'linestyle': '-', 'marker': '*', 'linewidth': 2.0, 'markersize': 10.0, 'label': 'omps_nm', 'fontsize': 14.4}
        if avg_window is None:
            # bug fixed (AttributeError: 'Rectangle' object has no property 'marker'). M.Li
            df[column].mean('y').plot.line(x = "time", ax=ax, color=plot_kwargs['color'],linestyle=plot_kwargs['linestyle'],\
            #df[column].mean('y').plot(ax=ax, color=plot_kwargs['color'],linestyle=plot_kwargs['linestyle'],\
                                           marker=plot_kwargs['marker'],linewidth=plot_kwargs['linewidth'],\
                                          markersize=plot_kwargs['markersize'],label=plot_kwargs['label'])
        else:
            # bug fixed (AttributeError: 'Rectangle' object has no property 'marker'). M.Li
            df[column].resample(time = avg_window).mean().mean('y').plot.line(x = "time", ax=ax,color=plot_kwargs['color'],\
            #df[column].resample(time = avg_window).mean().mean('y').plot(ax=ax,color=plot_kwargs['color'],\
                                                                              linestyle=plot_kwargs['linestyle'],\
                                           marker=plot_kwargs['marker'],linewidth=plot_kwargs['linewidth'],\
                                          markersize=plot_kwargs['markersize'],label=plot_kwargs['label'])
    
    # If plot has been created add to the current axes.
    else:
        # this means that an axis handle already exists and use it to plot the model output.
        if avg_window is None:
            # bug fixed. M.Li
            df[column].mean('y').plot.line(x = "time",ax=ax, color=plot_dict['color'],linestyle=plot_dict['linestyle'],\
            #df[column].mean('y').plot(ax=ax, color=plot_dict['color'],linestyle=plot_dict['linestyle'],\
                                           marker=plot_dict['marker'],linewidth=plot_dict['linewidth'],\
                                          markersize=plot_dict['markersize'],label=plot_dict['label'])
        else:
            # bug fixed. M.Li
            df[column].resample(time=avg_window).mean().mean('y').plot.line(x = "time",ax=ax, color=plot_dict['color'],\
            #df[column].resample(time=avg_window).mean().mean('y').plot(ax=ax, color=plot_dict['color'],\
                                                                            linestyle=plot_dict['linestyle'],\
                                           marker=plot_dict['marker'],linewidth=plot_dict['linewidth'],\
                                          markersize=plot_dict['markersize'],label=plot_dict['label'])   
    
    #Set parameters for all plots
    ax.set_ylabel(ylabel,fontweight='bold',**text_kwargs)
    ax.set_xlabel('time',fontweight='bold',**text_kwargs)
    ax.legend(frameon=False,fontsize=text_kwargs['fontsize']*0.8)
    ax.tick_params(axis='both',length=10.0,direction='inout')
    ax.tick_params(axis='both',which='minor',length=5.0,direction='out')
    ax.legend(frameon=False,fontsize=text_kwargs['fontsize']*0.8,
              bbox_to_anchor=(1.0, 0.9), loc='center left')
    if domain_type is not None and domain_name is not None:
        if domain_type == 'epa_region':
            ax.set_title('EPA Region ' + domain_name,fontweight='bold',**text_kwargs)
        else:
            ax.set_title(domain_name,fontweight='bold',**text_kwargs)
    return ax
    
def make_taylor(df,df_reg=None, column_o=None, label_o='Obs', column_m=None, label_m='Model', 
                dia=None, ylabel=None, ty_scale=1.5,
                domain_type=None, domain_name=None,
                plot_dict=None, fig_dict=None, text_dict=None,debug=False):
    """Creates taylor plot. Note sometimes model values are off the scale 
    on this plot. This will be fixed soon.
    
    Parameters
    ----------
    df : dataframe
        model/obs pair data to plot
    df_reg: not currently enabled. empty argument for symmetry with surfplots
        model/obs paired regulatory data to plot
    column_o : str
        Column label of observational variable to plot
    label_o : str
        Name of observational variable to use in plot legend
    column_m : str
        Column label of model variable to plot
    label_m : str
        Name of model variable to use in plot legend 
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
    nan_ind = ((~np.isnan(df[column_o].values))&(~np.isnan(df[column_m].values)))
    #First define items for all plots
    if debug is False:
        plt.ioff()
        
    #set default text size
    def_text = dict(fontsize=14.0)
    if text_dict is not None:
        text_kwargs = {**def_text, **text_dict}
    else:
        text_kwargs = def_text
    # set ylabel to column if not specified.
    if ylabel is None:
        ylabel = column_o
    #Then, if no plot has been created yet, create a plot and plot the first pair.
    if dia is None: 
        # create the figure
        if fig_dict is not None:
            f = plt.figure(**fig_dict)    
        else: 
            f = plt.figure(figsize=(12,10))    
        sns.set_style('ticks')
        # plot the line
        dia = td(df[column_o].std().values, scale=ty_scale, fig=f,
                               rect=111, label=label_o)
        plt.grid(linewidth=1, alpha=.5)
        cc = corrcoef(df[column_o].values[nan_ind].flatten(), df[column_m].values[nan_ind].flatten())[0, 1]
        dia.add_sample(df[column_m].std().values, cc, zorder=9, label=label_m, **plot_dict)
    # If plot has been created add to the current axes.
    else:
        # this means that an axis handle already exists and use it to plot another model
        cc = corrcoef(df[column_o].values[nan_ind].flatten(), df[column_m].values[nan_ind].flatten())[0, 1]
        dia.add_sample(df[column_m].std().values, cc, zorder=9, label=label_m, **plot_dict)
    #Set parameters for all plots
    contours = dia.add_contours(colors='0.5')
    # control the clabel format for very high values (e.g., NO2 columns), M.Li
    #plt.clabel(contours, inline=1, fontsize=text_kwargs['fontsize']*0.8)
    plt.clabel(contours, inline=1, fontsize=text_kwargs['fontsize']*0.8, fmt='(%1.1e)')

    plt.grid(alpha=.5)
    plt.legend(frameon=False,fontsize=text_kwargs['fontsize']*0.8,
               bbox_to_anchor=(0.75, 0.93), loc='center left')
    if domain_type is not None and domain_name is not None:
        if domain_type == 'epa_region':
            plt.title('EPA Region ' + domain_name,fontweight='bold',**text_kwargs)
        else:
            plt.title(domain_name,fontweight='bold',**text_kwargs)
    ax = plt.gca()
    ax.axis["left"].label.set_text('Standard Deviation: '+ylabel)
    ax.axis["top"].label.set_text('Correlation')
    ax.axis["left"].label.set_fontsize(text_kwargs['fontsize'])
    ax.axis["top"].label.set_fontsize(text_kwargs['fontsize'])
    ax.axis["left"].label.set_fontweight('bold')
    ax.axis["top"].label.set_fontweight('bold')
    ax.axis["top"].major_ticklabels.set_fontsize(text_kwargs['fontsize']*0.8)
    ax.axis["left"].major_ticklabels.set_fontsize(text_kwargs['fontsize']*0.8)
    ax.axis["right"].major_ticklabels.set_fontsize(text_kwargs['fontsize']*0.8)
    return dia

# expand multiboxplot capabilties to the satellite 
# multiboxplot and "diurnal" spatial overlays need this hour mask helper 
def _hour_window_mask(times, lon, hour_range, hour_basis="solar"):
    """Mask of samples whose hour-of-day falls inside ``hour_range``.

    Shared by the windowed spatial overlay so the
    window semantics stay identical. ``hour_range`` is end-exclusive and
    wraps midnight when start > end (e.g. [22, 4]). ``hour_basis``:
    "solar" = local solar time per cell (UTC + lon/15, so ``lon``
    broadcasts against time), "utc" = as-is, or a fixed UTC offset in
    hours (e.g. -6).

    Returns
    -------
    (mask, window_label) : (xr.DataArray of bool, str)
    """
    h0, h1 = (float(hour_range[0]), float(hour_range[1]))
    hh = times.dt.hour + times.dt.minute / 60.0
    b = str(hour_basis).lower()
    if b == "utc":
        lst = hh % 24
        basis_label = "UTC"
    elif b == "solar":
        lst = (hh + lon / 15.0) % 24
        basis_label = "LST"
    else:
        lst = (hh + float(hour_basis)) % 24
        basis_label = f"UTC{float(hour_basis):+g}h"
    if h0 <= h1:
        keep = (lst >= h0) & (lst < h1)
    else:  # window wraps midnight
        keep = (lst >= h0) | (lst < h1)
    return keep, f"{h0:g}-{h1:g} {basis_label}"
    
# driver calls splots.make_spatial_overlay and then resolves to sat plots. 
# sat plots does not have ucomp and vcomp. so silently fail it. 

def make_spatial_overlay(df, vmodel, column_o=None, label_o=None, column_m=None,
                      label_m=None, ylabel=None, vmin=None,
                      vmax=None, nlevels=None, proj=None, outname='plot',
                      u_comp=None, v_comp=None, wind_barb=False, wind_barb_step = 1,  wind_barb_kwargs=None,
                      domain_type=None, domain_name=None, fig_dict=None,
                      text_dict=None, debug=False, uxgrid=None, gridlines = False, reduction_dict = None):
        
    """Creates spatial overlay plot. 
    
    Parameters
    ----------
    df : dataframe
        model/obs pair data to plot
    vmodel: dataarray
        slice of model data to plot
    column_o : str
        Column label of observation variable to plot
    label_o : str
        Name of observation variable to use in plot title 
    column_m : str
        Column label of model variable to plot
    label_m : str
        Name of model variable to use in plot title
    ylabel : str
        Title of colorbar axis
    vmin : real number
        Min value to use on colorbar axis
    vmax : real number
        Max value to use on colorbar axis
    nlevels: integer
        Number of levels used in colorbar axis
    proj: cartopy projection
        cartopy projection to use in plot
    outname : str
        file location and name of plot (do not include .png)
    domain_type : str
        Domain type specified in input yaml file
    domain_name : str
        Domain name specified in input yaml file
    fig_dict : dictionary
        Dictionary containing information about figure
           - includes: cbar orientation and cbar kwargs
    text_dict : dictionary
        Dictionary containing information about text
    debug : boolean
        Whether to plot interactively (True) or not (False). Flag for 
        submitting jobs to supercomputer turn off interactive mode.
    gridlines: boolean 
        Draw lat lon lines with labels on each map 
    reduction_dict : dictionary, optional
        Controls how the per-granule paired time series is reduced to the
        single map (satellite Dataset path only). Set in the plot group's
        data_proc section of the YAML. Keys:
        - time_reduction : {"mean", "median"}, default "mean"
        - daily_first : bool, default False. First average granules within
          each day, then reduce across days, so every day carries equal
          weight regardless of how many valid scans/orbits it has.
        - common_mask : bool, default True. Only reduce over times where BOTH
          obs and model are valid, so all three panels (and the bias) are
          computed from identical samples.
        - min_obs : int, default 0. Mask cells with fewer than this many
          valid granule samples over the analysis window.
        - hour_range : [start, end], optional. Keep only samples whose
          hour-of-day falls in this window (end exclusive; wraps midnight
          if start > end, e.g. [22, 4]) before the time reduction. Use to
          build e.g. morning-commute maps from hourly TEMPO scans.
        - hour_basis : {"solar", "utc", <number>}, default "solar". How
          hour_range is evaluated: "solar" = local solar time per cell
          (UTC + lon/15), "utc" = as-is, or a fixed UTC offset in hours
          (e.g. -6). The window is annotated in the figure suptitle.
          
        
    Returns
    -------
    plot 
        spatial overlay plot
        
    """
    if debug is False:
        plt.ioff()
        
    def_map = dict(states=True,figsize=[15, 8])
    if fig_dict is not None:
        map_kwargs = {**def_map, **fig_dict}
    else:
        map_kwargs = def_map
  
    #set default text size
    def_text = dict(fontsize=20)
    if text_dict is not None:
        text_kwargs = {**def_text, **text_dict}
    else:
        text_kwargs = def_text
        
    # set ylabel to column if not specified.
    if ylabel is None:
        ylabel = column_o
    
    # #Take the mean for each siteid
    # df_mean=df.groupby(['siteid'],as_index=False).mean(numeric_only=True)

    # satellite paired data is xarray not a pandas df. need to guard against it
    _df_is_ds = hasattr(df, "data_vars")
    if _df_is_ds:
        df_mean = None
    else:
        df_mean = df.groupby(['siteid'], as_index=False).mean(numeric_only=True)
    
    #Take the mean over time for the model output
    vmodel_mean = vmodel[column_m].mean(dim='time').squeeze()
    
    #Determine the domain
    if domain_type == 'all' and domain_name == 'CONUS':
        latmin= 25.0
        lonmin=-130.0
        latmax= 50.0
        lonmax=-60.0
        title_add = domain_name + ': '
    elif domain_type == 'epa_region' and domain_name is not None:
        latmin,lonmin,latmax,lonmax,_ = get_epa_bounds(index=None,acronym=domain_name)
        title_add = 'EPA Region ' + domain_name + ': '
    else:
        # # float should work for both 1D coords for sat datasets and pandas
        # latmin = math.floor(float(df.latitude.min()))
        # lonmin = math.floor(float(df.longitude.min()))
        # latmax = math.ceil(float(df.latitude.max()))
        # lonmax = math.ceil(float(df.longitude.max()))
        # latmin= math.floor(min(df.latitude))
        # lonmin= math.floor(min(df.longitude))
        # latmax= math.ceil(max(df.latitude))
        # lonmax= math.ceil(max(df.longitude))
        # title_add = domain_name + ': '

        # zoom in where there is valid data rather than make a global map 
        try:
            _m = np.isfinite(df[column_o])
            for _d in list(_m.dims):
                if _d not in df["latitude"].dims:      # collapse time etc.
                    _m = _m.any(_d)
            _lat = df["latitude"].where(_m)
            _lon = df["longitude"].where(_m)
            latmin = math.floor(float(_lat.min())); lonmin = math.floor(float(_lon.min()))
            latmax = math.ceil(float(_lat.max())); lonmax = math.ceil(float(_lon.max()))
        except Exception:
            latmin = math.floor(float(df.latitude.min()))
            lonmin = math.floor(float(df.longitude.min()))
            latmax = math.ceil(float(df.latitude.max()))
            lonmax = math.ceil(float(df.longitude.max()))
        title_add = (domain_name + ': ') if domain_name else ''
    
    #Map the model output first.
    cbar_kwargs = dict(aspect=15,shrink=.8)
    
    #Add options that this could be included in the fig_kwargs in yaml file too.
    if 'extent' not in map_kwargs:
        map_kwargs['extent'] = [lonmin,lonmax,latmin,latmax] 
    if 'crs' not in map_kwargs:
        map_kwargs['crs'] = proj
    
    #With pcolormesh, a Warning shows because nearest interpolation may not work for non-monotonically increasing regions.
    #Because I do not want to pull in the edges of the lat lon for every model I switch to contourf.
    #First determine colorbar, so can use the same for both contourf and scatter
    
    if vmin is None and vmax is None:
        # obs source for color limits: df_mean[column_o] (station path) or
        # df[column_o] (sat Dataset path)
        _obs_for_limits = df_mean[column_o] if df_mean is not None else df[column_o]
        vmin = float(np.min((vmodel_mean.quantile(0.01), _obs_for_limits.quantile(0.01))))
        vmax = float(np.max((vmodel_mean.quantile(0.99), _obs_for_limits.quantile(0.99))))
        
        # vmin = np.min((vmodel_mean.quantile(0.01), df_mean[column_o].quantile(0.01)))
        # vmax = np.max((vmodel_mean.quantile(0.99), df_mean[column_o].quantile(0.99)))
        
    if nlevels is None:
        nlevels = 21
    
    clevel = np.linspace(vmin,vmax,nlevels)
    cmap = mpl.cm.get_cmap('Spectral_r',nlevels-1) 
    norm = mpl.colors.BoundaryNorm(clevel, ncolors=cmap.N, clip=False)
        
    # #I add extend='both' here because the colorbar is setup to plot the values outside the range
    # ax = vmodel_mean.monet.quick_contourf(cbar_kwargs=cbar_kwargs, figsize=map_kwargs['figsize'], map_kws=map_kwargs,
    #                             robust=True, norm=norm, cmap=cmap, levels=clevel, extend='both') 

    # Structured -> quick_contourf (2-D lat/lon mesh). Unstructured (CESM-SE
    # ncol/n_face) goes via monet.draw_map + uxarray PolyCollection so we
    if _df_is_ds:
        red = reduction_dict or {}
        how = str(red.get("time_reduction", "mean")).lower()
        if how not in ("mean", "median"):
            print(f"make_spatial_overlay: unknown time_reduction '{how}', using 'mean'.")
            how = "mean"
        daily_first = bool(red.get("daily_first", False))
        common_mask = bool(red.get("common_mask", True))
        min_obs = int(red.get("min_obs", 0))

        obs_da = df[column_o]
        mod_da = vmodel[column_m]

        _window_label = None # have a way to subset diurnal plots by a range 
        
        if "time" in obs_da.dims and "time" in mod_da.dims:
            # Optional diurnal window:
            # evaluated in UTC, at a fixed UTC offset,
            # or in local solar time (longitude-dependent, UTC + lon/15 
            # so the window means the same local hours across the whole map).
            
            hour_range = red.get("hour_range")
            if hour_range is not None:
                _keep, _window_label = _hour_window_mask(
                    obs_da["time"], df["longitude"], hour_range,
                    red.get("hour_basis", "solar"))
                
                obs_da = obs_da.where(_keep)
                mod_da = mod_da.where(_keep)

                if not bool(_keep.any()):
                    print(
                        f"make_spatial_overlay: no samples fall in hour_range "
                        f"{_window_label} (a polar orbiter like TROPOMI only "
                        "samples ~13:30 LST); map will be empty."
                    )
            if common_mask and set(obs_da.dims) == set(mod_da.dims):
                _valid = obs_da.notnull() & mod_da.notnull()
                obs_da = obs_da.where(_valid)
                mod_da = mod_da.where(_valid)
            # sample count per cell, before any daily compositing
            n_valid = obs_da.notnull().sum("time")
            if daily_first:
                obs_da = obs_da.resample(time="1D").mean()
                mod_da = mod_da.resample(time="1D").mean()
            obs_field = getattr(obs_da, how)("time")
            mod_field = getattr(mod_da, how)("time")
            if min_obs > 0:
                obs_field = obs_field.where(n_valid >= min_obs)
                mod_field = mod_field.where(n_valid >= min_obs)
        else:
            obs_field = obs_da
            mod_field = vmodel_mean

        obs_field = obs_field.squeeze()
        mod_field = mod_field.squeeze()
        
        diff_field = (mod_field - obs_field)

        _is_unstruct = uxgrid is not None or any(
            d in mod_field.dims for d in ("n_face", "ncol"))
        if _is_unstruct and uxgrid is None:
            _gf = vmodel.attrs.get("mio_scrip_file") or vmodel.attrs.get("mio_grid_file")
            if not _gf:
                raise ValueError(
                    "satplots.make_spatial_overlay: unstructured model but no "
                    "uxgrid passed and no mio_scrip_file/mio_grid_file attr."
                )
    
            uxgrid = ux.open_grid(_gf)
            
        if _is_unstruct:
            from melodies_monet.plots.uxarray_render import render_unstructured_field
            
            def _draw(ax_, fld, cm, nm):
                return render_unstructured_field(
                    ax_, fld, uxgrid, cmap=cm, norm=nm,
                    extent=map_kwargs["extent"], coast=True, borders=True,
                    states=map_kwargs.get("states", True), gridlines=gridlines, colorbar=False)
        else:
            import cartopy.feature as cfeature
            
            def _draw(ax_, fld, cm, nm):
                pm = ax_.pcolormesh(
                    np.asarray(fld["longitude"].values), np.asarray(fld["latitude"].values),
                    np.asarray(fld.values), cmap=cm, norm=nm,
                    transform=ccrs.PlateCarree(), shading="auto")
                ax_.coastlines(linewidth=0.5)
                ax_.add_feature(cfeature.BORDERS, linewidth=0.4)
                
                if map_kwargs.get("states", True):
                    ax_.add_feature(cfeature.STATES, linewidth=0.3)
                ax_.set_extent(map_kwargs["extent"], crs=ccrs.PlateCarree())
                
                if gridlines:
                    # same style as uxarray_render 
                    gl = ax_.gridlines(draw_labels=True, lw=1.0, color="black",
                                       alpha=0.5, linestyle=":")
                    gl.top_labels = False
                    gl.right_labels = False
                                    
                return pm
                
        _proj = proj if proj is not None else ccrs.PlateCarree()
        figsize = map_kwargs.get("figsize", [22, 6])
        
        # Colorbar layout
        # Orientation is a YAML option
        cbar_orientation = str(map_kwargs.pop("cbar_orientation", "vertical")).lower()
        _user_cbk = map_kwargs.pop("cbar_kwargs", None) or {}
        if cbar_orientation == "horizontal":
            cbk = dict(location="bottom", shrink=0.7, aspect=35, pad=0.04, extend="both")
        else:
            cbk = dict(location="right", shrink=0.75, aspect=25, pad=0.02, extend="both")
        cbk.update(_user_cbk)

        fig, axes = plt.subplots(1, 3, figsize=figsize,
                                 subplot_kw={"projection": _proj},
                                 constrained_layout=True)
        
        # obs (left) + model (middle): shared Spectral_r scale
        _draw(axes[0], obs_field, cmap, norm)
        axes[0].set_title(label_o, fontweight="bold", **text_kwargs)
        poly = _draw(axes[1], mod_field, cmap, norm)
        axes[1].set_title(label_m, fontweight="bold", **text_kwargs)
        cbar = fig.colorbar(poly, ax=axes[:2].tolist(), **cbk)
        cbar.set_label(ylabel, fontweight="bold", **text_kwargs)
        cbar.ax.tick_params(labelsize=text_kwargs["fontsize"] * 0.8)

        # bias (right): model - obs, diverging scale
        _dv = np.asarray(diff_field.values, dtype=float)
        _dv = _dv[np.isfinite(_dv)]
        _vd = float(np.nanmax(np.abs(np.percentile(_dv, [1, 99])))) if _dv.size else 1.0
        if not np.isfinite(_vd) or _vd == 0:
            _vd = 1.0
        _bn = nlevels if nlevels else 21
        _bcmap = mpl.cm.get_cmap("RdBu_r", _bn - 1)
        _bnorm = mpl.colors.BoundaryNorm(np.linspace(-_vd, _vd, _bn), ncolors=_bcmap.N, clip=False)
        bpoly = _draw(axes[2], diff_field, _bcmap, _bnorm)
        axes[2].set_title(label_m + " - " + label_o, fontweight="bold", **text_kwargs)

        bcbar = fig.colorbar(bpoly, ax=axes[2], ticks=np.linspace(-_vd, _vd, 5), **cbk)
        bcbar.set_label(ylabel, fontweight="bold", **text_kwargs)
        bcbar.ax.tick_params(labelsize=text_kwargs["fontsize"] * 0.8)
        _bax = bcbar.ax.xaxis if cbk.get("location") == "bottom" else bcbar.ax.yaxis
        _bax.set_major_formatter(mpl.ticker.FormatStrFormatter("%.2g"))

        _suptitle = (title_add or "").strip().rstrip(":").strip()
        if _window_label:
            _suptitle = f"{_suptitle} [{_window_label}]" if _suptitle else f"[{_window_label}]"
   
        if _suptitle:
            fig.suptitle(_suptitle, fontweight="bold", **text_kwargs)
        savefig(outname + ".png", loc=4, logo_height=100, bbox_inches="tight", dpi=150)

        if debug is False:
            plt.close(fig)  # long multi-group jobs otherwise accumulate open figures
                    
        return axes[1]
    
def calculate_boxplot(df, df_reg=None,column=None, label=None, plot_dict=None, comb_bx = None, label_bx = None):
    """Combines data into acceptable format for box-plot
    
    Parameters
    ----------
    df : dataframe
         model/obs pair data to plot
    df_reg: not currently enabled. empty argument for symmetry with surfplots
        model/obs paired regulatory data to plot
    column : str
        Column label of variable to plot
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
        #First define the colors for the observations.
        obs_dict = dict(color='gray', linestyle='-',marker='x', linewidth=1.2, markersize=6.)
        if plot_dict is not None:
            #Whatever is not defined in the yaml file is filled in with the obs_dict here.
            plot_kwargs = {**obs_dict, **plot_dict}
        else:
            plot_kwargs = obs_dict
    else:
        plot_kwargs = plot_dict
    #For all, a column to the dataframe and append the label info to the list.
    plot_kwargs['column'] = column
    plot_kwargs['label'] = label
    comb_bx[label] = df[column]
    label_bx.append(plot_kwargs)
    
    return comb_bx, label_bx
    
def make_boxplot(comb_bx, label_bx, ylabel = None, vmin = None, vmax = None, outname='plot',
                 domain_type=None, domain_name=None,
                 plot_dict=None, fig_dict=None,text_dict=None,debug=False,
                 set_stat_sig=False, gridlines=False):
    
    """Creates box-plot. 
    
    Parameters
    ----------
    comb_bx: dataframe
        dataframe containing information to create box-plot from 
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
    if debug is False:
        plt.ioff()
    #First define items for all plots
    #set default text size
    def_text = dict(fontsize=14)
    if text_dict is not None:
        text_kwargs = {**def_text, **text_dict}
    else:
        text_kwargs = def_text
    # set ylabel to column if not specified.
    if ylabel is None:
        ylabel = label_bx[0]
    
    #Fix the order and palate colors
    order_box = []
    pal = {}
    for i in range(len(label_bx)):
        order_box.append(label_bx[i]['label'])
        pal[label_bx[i]['label']] = label_bx[i]['color']
        
    #Make plot
    if fig_dict is not None:
        f,ax = plt.subplots(**fig_dict)    
    else: 
        f,ax = plt.subplots(figsize=(8,8))
    #Define characteristics of boxplot.
    boxprops = {'edgecolor': 'k', 'linewidth': 1.5}
    lineprops = {'color': 'k', 'linewidth': 1.5}
    boxplot_kwargs = {'boxprops': boxprops, 'medianprops': lineprops,
                  'whiskerprops': lineprops, 'capprops': lineprops,
                  'fliersize' : 2.0, 
                  'flierprops': dict(marker='*', 
                                     markerfacecolor='blue', 
                                     markeredgecolor='none',
                                     markersize = 6.0),
                  'width': 0.75, 'palette': pal,
                  'order': order_box,
                  'showmeans': True, 
                  'meanprops': {'marker': ".", 'markerfacecolor': 'black', 
                                'markeredgecolor': 'black',
                               'markersize': 20.0}}
    sns.set_style("whitegrid")
    sns.set_style("ticks")
    sns.boxplot(ax=ax,x="variable", y="value",data=pd.melt(comb_bx), **boxplot_kwargs)
    ax.set_xlabel('')
    ax.set_ylabel(ylabel,fontweight='bold',**text_kwargs)
    ax.tick_params(labelsize=text_kwargs['fontsize']*0.8)
    if domain_type is not None and domain_name is not None:
        if domain_type == 'epa_region':
            ax.set_title('EPA Region ' + domain_name,fontweight='bold',**text_kwargs)
        else:
            ax.set_title(domain_name,fontweight='bold',**text_kwargs)
    if vmin is not None and vmax is not None:
        ax.set_ylim(ymin = vmin, ymax = vmax)
    
    plt.tight_layout()
    savefig(outname + '.png',loc=4, logo_height=100, bbox_inches='tight', dpi=200)
    if debug is False:
        plt.close(plt.gcf())  # free the figure; long jobs accumulate otherwise
        
def make_spatial_bias_gridded(df, column_o=None, label_o=None, column_m=None, 
                      label_m=None, ylabel = None, vmin=None, vdiff=None,
                      vmax = None, nlevels = None, proj = None, outname = 'plot', 
                      domain_type=None, domain_name=None, fig_dict=None, 
                      text_dict=None,debug=False):
        
    """Creates difference plot for satellite and model data.
        For data in swath format, overplots all differences
        For data on regular grid, mean difference.
    """
    if debug is False:
        plt.ioff()
        
    def_map = dict(states=True,figsize=[15, 8])
    if fig_dict is not None:
        map_kwargs = {**def_map, **fig_dict}
    else:
        map_kwargs = def_map
  
    #set default text size
    def_text = dict(fontsize=20)
    if text_dict is not None:
        text_kwargs = {**def_text, **text_dict}
    else:
        text_kwargs = def_text
        
    # set ylabel to column if not specified.
    if ylabel is None:
        ylabel = column_o
    
    #Take the difference for the model output - the sat output

    diff_mod_min_obs = (df[column_m] - df[column_o]).squeeze()
    # #Take mean over time, 
    # if len(diff_mod_min_obs.dims) == 3:
    #     diff_mod_min_obs = diff_mod_min_obs.mean('time')

    # Reduce away time so the renderer (uxarray polygons or pcolormesh)
    # gets a single map regardless of timestep count.
    diff_mod_min_obs = df[column_m] - df[column_o]
    if "time" in diff_mod_min_obs.dims:
        diff_mod_min_obs = diff_mod_min_obs.mean("time")
    diff_mod_min_obs = diff_mod_min_obs.squeeze()
    
    #Determine the domain
    if domain_type == 'all' and domain_name == 'CONUS':
        latmin= 25.0
        lonmin=-130.0
        latmax= 50.0
        lonmax=-60.0
        title_add = domain_name + ': '
    elif domain_type == 'epa_region' and domain_name is not None:
        latmin,lonmin,latmax,lonmax,_ = get_epa_bounds(index=None,acronym=domain_name)
        title_add = 'EPA Region ' + domain_name + ': '
    else:
        # latmin= -90
        # lonmin= -180
        # latmax= 90
        # lonmax= 180
        title_add = (domain_name + ': ') if domain_name else ''
        # zoom to where the data actually is
        try:
            _m = np.isfinite(diff_mod_min_obs)
            _lon = df["longitude"].where(_m)
            _lat = df["latitude"].where(_m)
            lonmin = float(_lon.min()); lonmax = float(_lon.max())
            latmin = float(_lat.min()); latmax = float(_lat.max())
            _padx = max(0.5, 0.05 * (lonmax - lonmin))
            _pady = max(0.5, 0.05 * (latmax - latmin))
            lonmin -= _padx; lonmax += _padx
            latmin -= _pady; latmax += _pady
        except Exception:
            latmin, lonmin, latmax, lonmax = -90, -180, 90, 180        
    
    #Map the model output first.
    cbar_kwargs = dict(aspect=15,shrink=.8)
    
    #Add options that this could be included in the fig_kwargs in yaml file too.
    if 'extent' not in map_kwargs:
        map_kwargs['extent'] = [lonmin,lonmax,latmin,latmax] 
    if 'crs' not in map_kwargs:
        map_kwargs['crs'] = proj
    
    #First determine colorbar
    if vmin is None and vmax is None and vdiff is None:
        #vmin = vmodel_mean.quantile(0.01)
        vmax = np.max((np.abs(diff_mod_min_obs.quantile(0.99)),np.abs(diff_mod_min_obs.quantile(0.01))))
        vmin = -vmax
    if vdiff is not None:
        vmax = np.float64(vdiff)
        vmin = -np.float64(vdiff)
        
    if nlevels is None:
        nlevels = 21
    print(vmin,vmax)
    clevel = np.linspace(vmin,vmax,nlevels)
    cmap = mpl.cm.get_cmap('bwr',nlevels-1) 
    norm = mpl.colors.BoundaryNorm(clevel, ncolors=cmap.N, clip=False)
        
    #I add extend='both' here because the colorbar is setup to plot the values outside the range
    states = fig_dict.get('states', False)
    counties = fig_dict.get('counties', False)
    ax = monet.plots.mapgen.draw_map(crs=map_kwargs['crs'],extent=map_kwargs['extent'], states=states, counties=counties)
    # draw scatter plot of model and satellite differences
    markersize = fig_dict.get('markersize', 2)
    c = ax.axes.scatter(df.longitude,df.latitude,c=diff_mod_min_obs,cmap=cmap,s=markersize,norm=norm)
    plt.gcf().canvas.draw() 
    plt.tight_layout(pad=0)
    plt.title(title_add + label_m + ' - ' + label_o,fontweight='bold',**text_kwargs)
    ax.axes.set_extent(map_kwargs['extent'],crs=ccrs.PlateCarree())    
    
    #Uncomment these lines if you update above just to verify colorbars are identical.
    #Also specify plot above scatter = ax.axes.scatter etc.
    #cbar = ax.figure.get_axes()[1] 
    plt.colorbar(c,ax=ax,extend='both',**cbar_kwargs)
    
    #Update colorbar
    f = plt.gcf()
    
    model_ax = f.get_axes()[0]
    cax = f.get_axes()[1]
    
    #get the position of the plot axis and use this to rescale nicely the color bar to the height of the plot.
    position_m = model_ax.get_position()
    position_c = cax.get_position()
    cax.set_position([position_c.x0, position_m.y0, position_c.x1 - position_c.x0, (position_m.y1-position_m.y0)*1.1])
    cax.set_ylabel(r'$\Delta$'+ylabel,fontweight='bold',**text_kwargs)
    cax.tick_params(labelsize=text_kwargs['fontsize']*0.8,length=10.0,width=2.0,grid_linewidth=2.0)    
    
    #plt.tight_layout(pad=0)
    savefig(outname + '.png',loc=4, logo_height=100, bbox_inches='tight', dpi=150)
    if debug is False:
        plt.close(plt.gcf())  # free the figure; long jobs accumulate otherwise
    return ax    

def calculate_multi_boxplot(df, df_reg=None, region_name=None,
                            interval_list=None, interval_var=None,
                            interval_labels=None, column=None, label=None,
                            plot_dict=None, comb_bx=None, label_bx=None, hour_range=None, hour_basis="solar"):
    """Accumulate per-interval box statistics for one obs or model column.

    Satellite (xarray) counterpart of surfplots.calculate_multi_boxplot,
    with the same driver-facing signature. The x-axis bins come from
    ``interval_var`` + ``interval_list`` (bin edges):

    - ``interval_var: local_hour`` (or ``solar_hour``) -- local solar
      hour-of-day per cell (UTC + lon/15), e.g. edges [5, 11, 16, 20]
      with labels [morning, midday, evening]. The LST analog of the
      hour_range windowed spatial overlays.
    - ``interval_var: utc_hour`` (or ``hour``) -- hour-of-day in UTC.
    - any other name -- a variable in the paired dataset to bin on
      (e.g. the obs column itself, cloud fraction, ...), like the
      surface temperature-bin examples.

    ``df_reg`` and ``region_name`` are accepted for signature symmetry;
    EPA-region binning is not implemented on the satellite path.

    Returns
    -------
    (comb_bx, label_bx, region_bx)
        comb_bx : list of per-source dicts
        {label, color, bins: {label: {stats, n, blabel}}};
        label_bx : list of {column, label} (surface-compatible bookkeeping);
        region_bx : None (regions unsupported here).
    """
    if interval_list is None or len(interval_list) < 2:
        raise ValueError(
            "satellite multi_boxplot needs interval_list (bin edges) and "
            "interval_var; region-based multi_boxplot is not implemented "
            "for satellite pairs."
        )
    if region_name is not None:
        print("calculate_multi_boxplot: region binning is not implemented "
              "for satellite pairs; using interval_var bins.")
    if comb_bx is None:
        comb_bx = []
    if label_bx is None:
        label_bx = []

    edges = [float(e) for e in interval_list]
    nbin = len(edges) - 1
    if interval_labels is not None and len(interval_labels) != nbin:
        print("calculate_multi_boxplot: len(interval_labels) != number of "
              "bins; falling back to edge labels.")
        interval_labels = None
    bin_names = (list(interval_labels) if interval_labels is not None
                 else [f"{edges[i]:g}-{edges[i+1]:g}" for i in range(nbin)])

    da = df[column]

    # Optional hour-of-day window: clamp to a fixed overpass (e.g. TEMPO's
    # hourly data to the ~13:30 LST TROPOMI overpass)  
    # Masks the values to NaN outside the window; the finite-filter in the binning loop
    # then drops them. Time coord is preserved, so dayofweek/hour bins still work.
    if hour_range is not None:
        _keep, _win = _hour_window_mask(da["time"], df["longitude"], hour_range, hour_basis)
        da = da.where(_keep)
        print(f"calculate_multi_boxplot: clamped '{label}' to hour window {_win}.")
        
    iv = str(interval_var).lower()
    if iv in ("local_hour", "solar_hour"):
        hh = da["time"].dt.hour + da["time"].dt.minute / 60.0
        bv = (hh + df["longitude"] / 15.0) % 24
        unit = " LST"
    elif iv in ("utc_hour", "hour"):
        hh = da["time"].dt.hour + da["time"].dt.minute / 60.0
        bv = hh % 24
        unit = " UTC"
    elif interval_var in df:
        bv = df[interval_var]
        unit = ""
    elif iv in ("dayofweek", "weekday"):
        # Mon=0 .. Sun=6; one dim ('time'), broadcasts over the spatial dim in
        # da.where(mask). Well-defined for satellite (weekend NO2 effect).
        bv = da["time"].dt.dayofweek
        unit = ""
    else:
        raise KeyError(
            f"multi_boxplot interval_var '{interval_var}' is neither an "
            "hour keyword (local_hour/utc_hour) nor a variable in the "
            "paired dataset."
        )

    entry = {"label": label,
             "color": (plot_dict or {}).get("color", "0.7"),
             "bins": {}}
    for i in range(nbin):
        lo, hi = edges[i], edges[i + 1]
        mask = (bv >= lo) & ((bv <= hi) if i == nbin - 1 else (bv < hi))
        v = np.asarray(da.where(mask).values, dtype=float).ravel()
        v = v[np.isfinite(v)]
        blabel = f"{lo:g}-{hi:g}{unit}"
        if v.size == 0:
            print(f"calculate_multi_boxplot: no '{label}' samples in bin "
                  f"'{bin_names[i]}' [{blabel}]; box skipped.")
            continue
        p5, q1, med, q3, p95 = np.percentile(v, [5, 25, 50, 75, 95])
        entry["bins"][bin_names[i]] = dict(
            stats=dict(med=med, q1=q1, q3=q3, whislo=p5, whishi=p95,
                       mean=float(v.mean()), fliers=[]),
            n=int(v.size), blabel=blabel)
    comb_bx.append(entry)
    label_bx.append({"column": column, "label": label})
    return comb_bx, label_bx, None

def make_multi_boxplot(comb_bx, label_bx, region_bx=None, region_list=None,
                       interval_labels=None, model_name_list=None,
                       ylabel=None, xlabel=None, vmin=None, vmax=None,
                       outname="plot", domain_type=None, domain_name=None,
                       plot_dict=None, fig_dict=None, text_dict=None,
                       gridlines=False, debug=False):
    
    """Grouped boxplot: x-axis = interval bins, hue = data source.

    Satellite counterpart of surfplots.make_multi_boxplot (same
    driver-facing signature): within each bin the obs box and one box per
    model sit side by side. 
    
    """
    import matplotlib.patches as mpatches

    if debug is False:
        plt.ioff()
    text_kwargs = {**dict(fontsize=14), **(text_dict or {})}
    if ylabel is None:
        ylabel = label_bx[0]["column"] if label_bx else comb_bx[0]["label"]
    names = ([e["label"] for e in comb_bx] if model_name_list is None
             else list(model_name_list))

    # bin order: as given in interval_labels, else first-seen across sources
    if interval_labels is not None:
        order = [str(b) for b in interval_labels]
    else:
        order = []
        for e in comb_bx:
            for b in e["bins"]:
                if b not in order:
                    order.append(b)
    order = [b for b in order if any(b in e["bins"] for e in comb_bx)]
    if not order:
        print("make_multi_boxplot: no bin had data; no figure written.")
        return

    figsize = (fig_dict or {}).get("figsize", (10, 6))
    f, ax = plt.subplots(figsize=figsize)
    if gridlines:
        ax.grid(axis="y", alpha=0.5)
    nseries = len(comb_bx)
    group_w = 0.8
    box_w = group_w / nseries * 0.85
    meanprops = {"marker": ".", "markerfacecolor": "black",
                 "markeredgecolor": "black", "markersize": 14}
    for i, e in enumerate(comb_bx):
        stats, positions = [], []
        for j, b in enumerate(order):
            if b in e["bins"]:
                stats.append(e["bins"][b]["stats"])
                positions.append(
                    j - group_w / 2 + (i + 0.5) * group_w / nseries)
        if not stats:
            continue
        arts = ax.bxp(stats, positions=positions, widths=box_w,
                      showfliers=False, showmeans=True, patch_artist=True,
                      medianprops=dict(color="k", linewidth=1.5),
                      meanprops=meanprops)
        for bx_ in arts["boxes"]:
            bx_.set(facecolor=e["color"], alpha=0.7, edgecolor="k")

    def _btick(b):
        # annotate with the first source that has this bin (usually obs)
        for e in comb_bx:
            if b in e["bins"]:
                info = e["bins"][b]
                return f"{b}\n[{info['blabel']}]\nn={info['n']:,}"
        return b

    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([_btick(b) for b in order],
                       fontsize=text_kwargs["fontsize"] * 0.75)
    ax.set_xlim(-0.6, len(order) - 0.4)
    if xlabel:
        ax.set_xlabel(xlabel, fontweight="bold", **text_kwargs)
    ax.set_ylabel(ylabel, fontweight="bold", **text_kwargs)
    ax.tick_params(axis="y", labelsize=text_kwargs["fontsize"] * 0.8)
    if vmin is not None and vmax is not None:
        ax.set_ylim(vmin, vmax)
    ax.legend(handles=[mpatches.Patch(facecolor=e["color"], alpha=0.7,
                                      edgecolor="k", label=nm)
                       for e, nm in zip(comb_bx, names)],
              frameon=False, fontsize=text_kwargs["fontsize"] * 0.8)
    if domain_type is not None and domain_name is not None:
        if domain_type == "epa_region":
            ax.set_title("EPA Region " + str(domain_name),
                         fontweight="bold", **text_kwargs)
        else:
            ax.set_title(str(domain_name), fontweight="bold", **text_kwargs)

    plt.tight_layout()
    savefig(outname + ".png", loc=4, logo_height=100, bbox_inches="tight",
            dpi=200)
    if debug is False:
        plt.close(f)

def _swath_extent(lon, lat, pad=0.05):
    """[W, E, S, N] bounding box of the pixel cloud, padded a touch."""
    lo, hi = float(np.nanmin(lon)), float(np.nanmax(lon))
    la, lb = float(np.nanmin(lat)), float(np.nanmax(lat))
    dx = max((hi - lo) * pad, 0.02)
    dy = max((lb - la) * pad, 0.02)
    return [lo - dx, hi + dx, la - dy, lb + dy]

def _swath_binned_mean(lon, lat, val, lon_edges, lat_edges):
    """Mean of ``val`` per (lon,lat) cell to (nlat, nlon) array, NaN where empty."""
    ok = np.isfinite(lon) & np.isfinite(lat) & np.isfinite(val)
    s, _, _ = np.histogram2d(lon[ok], lat[ok], bins=[lon_edges, lat_edges],
                             weights=val[ok])
    c, _, _ = np.histogram2d(lon[ok], lat[ok], bins=[lon_edges, lat_edges])
    with np.errstate(invalid="ignore"):
        mean = s / c
    return np.where(c > 0, mean, np.nan).T 
    
def plot_swath_scatter(ds, model_var, obs_var, unc_var=None,
                       label_m="model", label_o="obs", ylabel=None,
                       outname="swath_scatter", extent=None, proj=None,
                       vmin=None, vmax=None, markersize=6, text_dict=None,
                       states=True, gridlines=True, render="auto", share_scale=True,
                       bin_deg=None, common_mask=True, debug=False): # create a lot of yaml customizations
    
    """Native-TEMPO pixel scatter: obs | model | bias(model-obs), 3 panels.

    Plots every swath pixel at its own longitude/latitude with ``ax.scatter``
    (no gridding). Operates on the 'swath' pair vector (dims ``obs``, coords
    ``longitude``/``latitude``). ``unc_var`` is unused here but accepted so
    callers can pass the pair's uncertainty name uniformly.

    Parameters mirror make_spatial_overlay where sensible. ``vmin``/``vmax``
    bound the obs+model color scale; if None they are the 2/98th percentiles
    of the obs values.
    """
    import cartopy.feature as cfeature

    text_kwargs = dict(text_dict) if text_dict else {"fontsize": 14}
    text_kwargs.setdefault("fontsize", 14)

    lon = np.asarray(ds["longitude"].values, dtype=float).ravel()
    lat = np.asarray(ds["latitude"].values, dtype=float).ravel()
    o = np.asarray(ds[obs_var].values, dtype=float).ravel()
    m = np.asarray(ds[model_var].values, dtype=float).ravel()

    # Common mask: the model is defined wherever the mesh covers, but the obs
    # is QA/cloud-filtered
    # eventually shouldnt this default to using the paired df which drops rows where there are nans in either 
    # obs / mod? 
    if common_mask:
        _both = np.isfinite(o) & np.isfinite(m)
        o = np.where(_both, o, np.nan)
        m = np.where(_both, m, np.nan)
        
    d = m - o
    n = int(np.isfinite(o).sum())

    if extent is None:
        extent = _swath_extent(lon, lat)
    _proj = proj if proj is not None else ccrs.PlateCarree()

    # true pixel footprints require the corner bounds carried on the vector
    _has_bounds = ("longitude_bounds" in ds and "latitude_bounds" in ds)
    
    mode = render
    if mode == "auto":
        # footprints (no resampling) are the correct native rendering when the
        # corner bounds are present; else fall back to scatter (small N) / bin.
        mode = "footprint" if _has_bounds else ("binned" if n > 5_000 else "scatter")
    if mode == "footprint" and not _has_bounds:
        print("plot_swath_scatter: render='footprint' but the swath has no "
              "longitude_bounds/latitude_bounds (re-pair to carry them); "
              "falling back to scatter.", flush=True)
        mode = "binned" if n > 5_000 else "scatter"
    if mode == "footprint" and n > 800_000:
        print(f"plot_swath_scatter: {n} footprint polygons is heavy and will "
              "overlap across overpasses -- consider time-subsetting to one "
              "overpass, or use the conservative obsgrid for aggregate maps.",
              flush=True)

    # Adaptive bin size when unset: target ~4 pixels/cell so the field stays
    # filled for sparse/coarse instruments (TROPOMI ~0.04-0.1 deg) 
    if bin_deg is None:
        _area = max((extent[1] - extent[0]) * (extent[3] - extent[2]), 1e-6)
        bin_deg = float(np.clip(np.sqrt(_area / max(n / 4.0, 1.0)), 0.02, 0.2))
        
    _finite = o[np.isfinite(o)]
    if vmin is None or vmax is None:
        _lo, _hi = (np.percentile(_finite, [2, 98]) if _finite.size else (0.0, 1.0))
        
        vmin = _lo if vmin is None else vmin
        vmax = _hi if vmax is None else vmax
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    cmap = mpl.cm.get_cmap("Spectral_r")
    _df = d[np.isfinite(d)]
    _vd = float(np.nanmax(np.abs(np.percentile(_df, [1, 99])))) if _df.size else 1.0
    if not np.isfinite(_vd) or _vd == 0:
        _vd = 1.0
    bnorm = mpl.colors.Normalize(vmin=-_vd, vmax=_vd)
    bcmap = mpl.cm.get_cmap("RdBu_r")

    if mode == "footprint":
        # draw each pixel as its true quadrilateral from the corner bounds --
        # no resampling
        from matplotlib.collections import PolyCollection
        _lonb = np.asarray(ds["longitude_bounds"].values, dtype=float)   # (obs, 4)
        _latb = np.asarray(ds["latitude_bounds"].values, dtype=float)
        _verts = np.stack([_lonb, _latb], axis=-1)                       # (obs, 4, 2)
        _ok = np.isfinite(_verts).all(axis=(1, 2))                       # complete quads

        def _paint(ax_, V, cm, nm):
            pc = PolyCollection(list(_verts[_ok]), array=np.asarray(V, float)[_ok],
                                cmap=cm, norm=nm, edgecolors="none", linewidths=0,
                                transform=ccrs.PlateCarree())
            pc.set_rasterized(True)
            ax_.add_collection(pc)
            return pc
        Fo, Fm, Fd = o, m, d
    elif mode == "binned":
        lon_e = np.arange(extent[0], extent[1] + bin_deg, bin_deg)
        lat_e = np.arange(extent[2], extent[3] + bin_deg, bin_deg)
        Fo = _swath_binned_mean(lon, lat, o, lon_e, lat_e)
        Fm = _swath_binned_mean(lon, lat, m, lon_e, lat_e)
        Fd = _swath_binned_mean(lon, lat, d, lon_e, lat_e)

        def _paint(ax_, F, cm, nm):
            return ax_.pcolormesh(lon_e, lat_e, F, cmap=cm, norm=nm,
                                  transform=ccrs.PlateCarree(), shading="flat")
    else:
        _s = markersize if n < 20_000 else (2 if n < 200_000 else 0.5)

        def _paint(ax_, V, cm, nm):
            return ax_.scatter(lon, lat, c=V, s=_s, marker="s", cmap=cm, norm=nm,
                               transform=ccrs.PlateCarree(), linewidths=0,
                               rasterized=True)
        Fo, Fm, Fd = o, m, d
        
    figsize = [22, 6]
    fig, axes = plt.subplots(1, 3, figsize=figsize,
                             subplot_kw={"projection": _proj},
                             constrained_layout=True)

    def _base(ax_):
        ax_.coastlines(linewidth=0.5)
        ax_.add_feature(cfeature.BORDERS, linewidth=0.4)
        
        if states:
            ax_.add_feature(cfeature.STATES, linewidth=0.3)
        if gridlines:
            gl = ax_.gridlines(draw_labels=True, lw=0.5, color="gray",
                               alpha=0.4, linestyle=":")
            gl.top_labels = False
            gl.right_labels = False
            
        ax_.set_extent(extent, crs=ccrs.PlateCarree())
    
    # model color scale: shared with obs (default)   
    if share_scale:
        norm_m = norm
    else:
        _fm = m[np.isfinite(m)]
        _mlo, _mhi = (np.percentile(_fm, [2, 98]) if _fm.size else (vmin, vmax))
        norm_m = mpl.colors.Normalize(vmin=_mlo, vmax=_mhi)
        
    im = _paint(axes[0], Fo, cmap, norm)
    
    _base(axes[0]); axes[0].set_title(label_o, fontweight="bold", **text_kwargs)
    # _paint(axes[1], Fm, cmap, norm)

    im_m = _paint(axes[1], Fm, cmap, norm_m)
    _base(axes[1]); axes[1].set_title(label_m, fontweight="bold", **text_kwargs)
    # cbar = fig.colorbar(im, ax=axes[:2].tolist(), location="right",
    #                     shrink=0.75, aspect=25, pad=0.02, extend="both")
    # cbar.set_label(ylabel, fontweight="bold", **text_kwargs)
    # cbar.ax.tick_params(labelsize=text_kwargs["fontsize"] * 0.8)

    if share_scale:
        cbar = fig.colorbar(im, ax=axes[:2].tolist(), location="right",
                            shrink=0.75, aspect=25, pad=0.02, extend="both")
        cbar.set_label(ylabel, fontweight="bold", **text_kwargs)
        cbar.ax.tick_params(labelsize=text_kwargs["fontsize"] * 0.8)
    else:
        # obs and model each get their own colorbar (different ranges)
        for _ax, _im in ((axes[0], im), (axes[1], im_m)):
            _cb = fig.colorbar(_im, ax=_ax, location="right",
                               shrink=0.75, aspect=25, pad=0.02, extend="both")
            _cb.set_label(ylabel, fontweight="bold", **text_kwargs)
            _cb.ax.tick_params(labelsize=text_kwargs["fontsize"] * 0.8)
            
    bim = _paint(axes[2], Fd, bcmap, bnorm)
    _base(axes[2])
    axes[2].set_title(label_m + " - " + label_o, fontweight="bold", **text_kwargs)
    bcbar = fig.colorbar(bim, ax=axes[2], location="right",
                         shrink=0.75, aspect=25, pad=0.02, extend="both")
    bcbar.set_label(ylabel, fontweight="bold", **text_kwargs)
    bcbar.ax.tick_params(labelsize=text_kwargs["fontsize"] * 0.8)

    _how = {"binned": f"{bin_deg:g}° binned", "footprint": "pixel footprints"}.get(
        mode, "native pixels")
    fig.suptitle(f"native swath ({_how}, n={n})", fontweight="bold", **text_kwargs)
    savefig(outname + ".png", loc=4, logo_height=100, bbox_inches="tight", dpi=150)
    if debug is False:
        plt.close(fig)

def plot_swath_obs_vs_model(ds, model_var, obs_var, unc_var=None,
                            label_m="model", label_o="obs", ylabel=None,
                            outname="swath_obs_vs_model", text_dict=None,
                            color_by_unc=True, xylim=None, debug=False):
    """Pixel-level obs-vs-model scatter for a native swath pair.

    The purest validation in the pipeline: every swath pixel already carries
    the observation and the AK-applied model colocated at the native
    footprint and overpass time, so this is a straight ``x = obs, y = model``
    scatter over the ``obs`` pixel dimension -- NO regridding, NO time
    averaging, NO colocation error. Draws the 1:1 line, a least-squares fit,
    and an N/r/slope/MB/RMSE box (computed on ALL finite pixels). If the pair
    carries per-pixel uncertainty (``unc_var``) the points are colored by it
    -- something the gridded products cannot show.

    Pass a single-overpass subset (e.g. from plot_one_overpass's TIME=/MATCH=
    selection) for a per-overpass scatter, or the full file for the aggregate.
    """
    tk = {"fontsize": 14, **(text_dict or {})}
    x = np.asarray(ds[obs_var].values, dtype=float).ravel()
    y = np.asarray(ds[model_var].values, dtype=float).ravel()
    c = None
    if color_by_unc and unc_var and unc_var in ds.variables:
        c = np.asarray(ds[unc_var].values, dtype=float).ravel()

    good = np.isfinite(x) & np.isfinite(y)
    if c is not None:
        c = c[good]
    x, y = x[good], y[good]
    if x.size < 3:
        print(f"plot_swath_obs_vs_model: only {x.size} valid pixels; skipping.",
              flush=True)
        return

    n = int(x.size)
    mb = float(np.mean(y - x))
    rmse = float(np.sqrt(np.mean((y - x) ** 2)))
    r = float(np.corrcoef(x, y)[0, 1])
    slope, intercept = (float(v) for v in np.polyfit(x, y, 1))
    
    if xylim is not None:
        lo, hi = float(xylim[0]), float(xylim[1])
    else:
        lo = float(min(x.min(), y.min()))
        hi = float(max(x.max(), y.max()))
        pad = 0.05 * (hi - lo or 1.0)
        lo, hi = lo - pad, hi + pad

    fig, ax = plt.subplots(figsize=(7.5, 7))
    ax.plot([lo, hi], [lo, hi], color="0.4", lw=1.2, zorder=1, label="1:1")
    sc = ax.scatter(x, y, c=(c if c is not None else "steelblue"),
                    s=8, alpha=0.5, edgecolors="none",
                    cmap="viridis" if c is not None else None, zorder=2)
    if c is not None:
        cb = fig.colorbar(sc, ax=ax, shrink=0.85)
        cb.set_label(f"{unc_var}", fontsize=tk["fontsize"] * 0.8)
    xs = np.linspace(lo, hi, 50)
    ax.plot(xs, slope * xs + intercept, "k--", lw=1.8, zorder=3,
            label=f"fit: slope={slope:.2f}")
    ax.text(0.03, 0.97,
            (f"N = {n}\nr = {r:.2f}\nslope = {slope:.2f}\n"
             f"MB = {mb:.2e}\nRMSE = {rmse:.2e}"),
            transform=ax.transAxes, va="top", ha="left",
            fontsize=tk["fontsize"] * 0.8,
            bbox=dict(boxstyle="round", fc="white", ec="0.6", alpha=0.85))
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect("equal", adjustable="box")
    yl = ylabel or "molecules/cm2"
    ax.set_xlabel(f"{label_o} ({yl})", fontweight="bold", fontsize=tk["fontsize"])
    ax.set_ylabel(f"{label_m} ({yl})", fontweight="bold", fontsize=tk["fontsize"])
    ax.legend(loc="lower right", fontsize=tk["fontsize"] * 0.75)
    savefig(outname + ".png", loc=4, logo_height=100, bbox_inches="tight", dpi=150)
    if debug is False:
        plt.close(fig)
    print(f"plot_swath_obs_vs_model: wrote {outname}.png (N={n}, r={r:.2f}, "
          f"slope={slope:.2f})", flush=True)

def plot_swath_oversampling(ds, bin_deg=0.02, outname="swath_oversampling",
                            extent=None, proj=None, text_dict=None,
                            states=True, debug=False):
    """Pixel-count density map: where TEMPO oversamples.

    2-D histogram of TEMPO pixel centers into ``bin_deg`` cells over the
    analysis window lead to repeat coverage / oversampling structure. Operates on
    the 'swath' pair vector (dims ``obs``, coords ``longitude``/``latitude``).
    """
    import cartopy.feature as cfeature

    text_kwargs = dict(text_dict) if text_dict else {"fontsize": 14}
    text_kwargs.setdefault("fontsize", 14)

    lon = np.asarray(ds["longitude"].values, dtype=float).ravel()
    lat = np.asarray(ds["latitude"].values, dtype=float).ravel()
    ok = np.isfinite(lon) & np.isfinite(lat)
    lon, lat = lon[ok], lat[ok]

    if extent is None:
        extent = _swath_extent(lon, lat)
    _proj = proj if proj is not None else ccrs.PlateCarree()

    lon_edges = np.arange(extent[0], extent[1] + bin_deg, bin_deg)
    lat_edges = np.arange(extent[2], extent[3] + bin_deg, bin_deg)
    counts, _, _ = np.histogram2d(lon, lat, bins=[lon_edges, lat_edges])
    counts = counts.T                      # (lat, lon) for pcolormesh
    counts = np.where(counts > 0, counts, np.nan)

    fig, ax = plt.subplots(figsize=[10, 8],
                           subplot_kw={"projection": _proj},
                           constrained_layout=True)
    cmap = mpl.cm.get_cmap("viridis")
    pm = ax.pcolormesh(lon_edges, lat_edges, counts, cmap=cmap,
                       transform=ccrs.PlateCarree(), shading="flat")
    ax.coastlines(linewidth=0.5)
    ax.add_feature(cfeature.BORDERS, linewidth=0.4)
    if states:
        ax.add_feature(cfeature.STATES, linewidth=0.3)
    ax.set_extent(extent, crs=ccrs.PlateCarree())
    cbar = fig.colorbar(pm, ax=ax, shrink=0.8, aspect=25, pad=0.02, extend="max")
    cbar.set_label(f"pixel count per {bin_deg:g}° cell",
                   fontweight="bold", **text_kwargs)
    cbar.ax.tick_params(labelsize=text_kwargs["fontsize"] * 0.8)
    ax.set_title(f"TEMPO sampling density (n={lon.size} pixels)",
                 fontweight="bold", **text_kwargs)
    savefig(outname + ".png", loc=4, logo_height=100, bbox_inches="tight", dpi=150)
    if debug is False:
        plt.close(fig)
    return ax
