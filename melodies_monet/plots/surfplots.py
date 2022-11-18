# Copyright (C) 2022 National Center for Atmospheric Research and National Oceanic and Atmospheric Administration
# SPDX-License-Identifier: Apache-2.0
#
#Code to create plots for surface observations

import os
import monetio as mio
import monet as monet
import seaborn as sns
from monet.util.tools import calc_8hr_rolling_max, calc_24hr_ave
import xarray as xr
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
from ..plots import savefig

# from util import write_ncf

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
    elif f.model.lower() == 'rrfs':
        proj = ccrs.LambertConformal(
            central_longitude=f.obj.cen_lon, central_latitude=f.obj.cen_lat)
    elif f.model.lower() in ['cesm_fv','cesm_se']:
        proj = ccrs.PlateCarree()
    elif f.model.lower() == 'random':
        proj = ccrs.PlateCarree()
    else: #Let's change this tomorrow to just plot as lambert conformal if nothing provided.
        raise NotImplementedError('Projection not defined for new model. Please add to surfplots.py')
    return proj

def make_spatial_bias(df, column_o=None, label_o=None, column_m=None, 
                      label_m=None, ylabel = None, vdiff=None,
                      outname = 'plot', 
                      domain_type=None, domain_name=None, fig_dict=None, 
                      text_dict=None,debug=False):
        
    """Creates surface spatial bias plot. 
    
    Parameters
    ----------
    df : dataframe
        model/obs pair data to plot
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
    vdiff : real number
        Min and max value to use on colorbar axis
    outname : str
        file location and name of plot (do not include .png)
    domain_type : str
        Domain type specified in input yaml file
    domain_name : str
        Domain name specified in input yaml file
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
        surface bias plot
        
    """
    if debug == False:
        plt.ioff()
        
    def_map = dict(states=True,figsize=[10, 5])
    if fig_dict is not None:
        map_kwargs = {**def_map, **fig_dict}
    else:
        map_kwargs = def_map
        
    #If not specified use the PlateCarree projection
    if 'crs' not in map_kwargs:
        map_kwargs['crs'] = ccrs.PlateCarree()
  
    #set default text size
    def_text = dict(fontsize=20)
    if text_dict is not None:
        text_kwargs = {**def_text, **text_dict}
    else:
        text_kwargs = def_text
        
    # set ylabel to column if not specified.
    if ylabel is None:
        ylabel = column_o
    
    #Take the mean for each siteid
    df_mean=df.groupby(['siteid'],as_index=False).mean()
       
    #Specify val_max = vdiff. the sp_scatter_bias plot in MONET only uses the val_max value
    #and then uses -1*val_max value for the minimum.
    ax = monet.plots.sp_scatter_bias(
        df_mean, col1=column_o, col2=column_m, map_kwargs=map_kwargs,val_max=vdiff,
        cmap=new_color_map(), edgecolor='k',linewidth=.8)
    
    if domain_type == 'all':
        latmin= 25.0
        lonmin=-130.0
        latmax= 50.0
        lonmax=-60.0
        plt.title(domain_name + ': ' + label_m + ' - ' + label_o,fontweight='bold',**text_kwargs)
    elif domain_type == 'epa_region' and domain_name is not None:
        latmin,lonmin,latmax,lonmax,acro = get_epa_bounds(index=None,acronym=domain_name)
        plt.title('EPA Region ' + domain_name + ': ' + label_m + ' - ' + label_o,fontweight='bold',**text_kwargs)
    else:
        latmin= math.floor(min(df.latitude))
        lonmin= math.floor(min(df.longitude))
        latmax= math.ceil(max(df.latitude))
        lonmax= math.ceil(max(df.longitude))
        plt.title(domain_name + ': ' + label_m + ' - ' + label_o,fontweight='bold',**text_kwargs)
    
    if 'extent' not in map_kwargs:
        map_kwargs['extent'] = [lonmin,lonmax,latmin,latmax]  
    ax.axes.set_extent(map_kwargs['extent'],crs=ccrs.PlateCarree())
    
    #Update colorbar
    f = plt.gcf()
    model_ax = f.get_axes()[0]
    cax = f.get_axes()[1]
    #get the position of the plot axis and use this to rescale nicely the color bar to the height of the plot.
    position_m = model_ax.get_position()
    position_c = cax.get_position()
    cax.set_position([position_c.x0, position_m.y0, position_c.x1 - position_c.x0, (position_m.y1-position_m.y0)*1.1])
    cax.set_ylabel(ylabel,fontweight='bold',**text_kwargs)
    cax.tick_params(labelsize=text_kwargs['fontsize']*0.8,length=10.0,width=2.0,grid_linewidth=2.0)    
    
    #plt.tight_layout(pad=0)
    savefig(outname + '.png', loc=4, logo_height=120)
    
def make_timeseries(df, column=None, label=None, ax=None, avg_window=None, ylabel=None,
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
    label : str
        Name of variable to use in plot legend 
    ax : ax
        matplotlib ax from previous occurrence so can overlay obs and model 
        results on the same plot
    avg_window : rule 
        Pandas resampling rule (e.g., 'H', 'D')
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
    if debug == False:
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
        if avg_window is None:
            ax = df[column].plot(ax=ax, **plot_kwargs)
        else:
            ax = df[column].resample(avg_window).mean().plot(ax=ax, legend=True, **plot_kwargs)
    
    # If plot has been created add to the current axes.
    else:
        # this means that an axis handle already exists and use it to plot the model output.
        if avg_window is None:
            ax = df[column].plot(ax=ax, legend=True, **plot_dict)
        else:
            ax = df[column].resample(avg_window).mean().plot(ax=ax, legend=True, **plot_dict)    
    
    #Set parameters for all plots
    ax.set_ylabel(ylabel,fontweight='bold',**text_kwargs)
    ax.set_xlabel(df.index.name,fontweight='bold',**text_kwargs)
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
    
def make_taylor(df, column_o=None, label_o='Obs', column_m=None, label_m='Model', 
                dia=None, ylabel=None, ty_scale=1.5,
                domain_type=None, domain_name=None,
                plot_dict=None, fig_dict=None, text_dict=None,debug=False):
    """Creates taylor plot. Note sometimes model values are off the scale 
    on this plot. This will be fixed soon.
    
    Parameters
    ----------
    df : dataframe
        model/obs pair data to plot
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
    #First define items for all plots
    if debug == False:
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
        dia = td(df[column_o].std(), scale=ty_scale, fig=f,
                               rect=111, label=label_o)
        plt.grid(linewidth=1, alpha=.5)
        cc = corrcoef(df[column_o].values, df[column_m].values)[0, 1]
        dia.add_sample(df[column_m].std(), cc, zorder=9, label=label_m, **plot_dict)
    # If plot has been created add to the current axes.
    else:
        # this means that an axis handle already exists and use it to plot another model
        cc = corrcoef(df[column_o].values, df[column_m].values)[0, 1]
        dia.add_sample(df[column_m].std(), cc, zorder=9, label=label_m, **plot_dict)
    #Set parameters for all plots
    contours = dia.add_contours(colors='0.5')
    plt.clabel(contours, inline=1, fontsize=text_kwargs['fontsize']*0.8)
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

def make_spatial_overlay(df, vmodel, column_o=None, label_o=None, column_m=None, 
                      label_m=None, ylabel = None, vmin=None,
                      vmax = None, nlevels = None, proj = None, outname = 'plot', 
                      domain_type=None, domain_name=None, fig_dict=None, 
                      text_dict=None,debug=False):
        
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
    text_dict : dictionary
        Dictionary containing information about text
    debug : boolean
        Whether to plot interactively (True) or not (False). Flag for 
        submitting jobs to supercomputer turn off interactive mode.
        
    Returns
    -------
    plot 
        spatial overlay plot
        
    """
    if debug == False:
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
    
    #Take the mean for each siteid
    df_mean=df.groupby(['siteid'],as_index=False).mean()
    
    #Take the mean over time for the model output
    vmodel_mean = vmodel[column_m].mean(dim='time').squeeze()
    
    #Determine the domain
    if domain_type == 'all':
        latmin= 25.0
        lonmin=-130.0
        latmax= 50.0
        lonmax=-60.0
        title_add = domain_name + ': '
    elif domain_type == 'epa_region' and domain_name is not None:
        latmin,lonmin,latmax,lonmax,acro = get_epa_bounds(index=None,acronym=domain_name)
        title_add = 'EPA Region ' + domain_name + ': '
    else:
        latmin= math.floor(min(df.latitude))
        lonmin= math.floor(min(df.longitude))
        latmax= math.ceil(max(df.latitude))
        lonmax= math.ceil(max(df.longitude))
        title_add = domain_name + ': '
    
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
    if vmin == None and vmax == None:
        vmin = np.min((vmodel_mean.quantile(0.01), df_mean[column_o].quantile(0.01)))
        vmax = np.max((vmodel_mean.quantile(0.99), df_mean[column_o].quantile(0.99)))
        
    if nlevels == None:
        nlevels = 21
    
    clevel = np.linspace(vmin,vmax,nlevels)
    cmap = mpl.cm.get_cmap('Spectral_r',nlevels-1) 
    norm = mpl.colors.BoundaryNorm(clevel, ncolors=cmap.N, clip=False)
        
    # For unstructured grid, we need a more advanced plotting code
    # Call an external funtion (Plot_2D)
    if vmodel.attrs.get('mio_has_unstructured_grid',False):
        from .Plot_2D import Plot_2D
        
        fig = plt.figure( figsize=fig_dict['figsize'] )
        ax = fig.add_subplot(1,1,1,projection=proj)
        
        p2d = Plot_2D( vmodel_mean, scrip_file=vmodel.mio_scrip_file, cmap=cmap, #colorticks=clevel, colorlabels=clevel,
                       cmin=vmin, cmax=vmax, lon_range=[lonmin,lonmax], lat_range=[latmin,latmax],
                       ax=ax, state=fig_dict['states'] )
    else:
        #I add extend='both' here because the colorbar is setup to plot the values outside the range
        ax = vmodel_mean.monet.quick_contourf(cbar_kwargs=cbar_kwargs, figsize=map_kwargs['figsize'], map_kws=map_kwargs,
                                    robust=True, norm=norm, cmap=cmap, levels=clevel, extend='both') 
    
    
    plt.gcf().canvas.draw() 
    plt.tight_layout(pad=0)
    plt.title(title_add + label_o + ' overlaid on ' + label_m,fontweight='bold',**text_kwargs)
     
    ax.axes.scatter(df_mean.longitude.values, df_mean.latitude.values,s=30,c=df_mean[column_o], 
                    transform=ccrs.PlateCarree(), edgecolor='b', linewidth=.50, norm=norm, 
                    cmap=cmap)
    ax.axes.set_extent(map_kwargs['extent'],crs=ccrs.PlateCarree())    
    
    #Uncomment these lines if you update above just to verify colorbars are identical.
    #Also specify plot above scatter = ax.axes.scatter etc.
    #cbar = ax.figure.get_axes()[1] 
    #plt.colorbar(scatter,ax=ax)
    
    #Update colorbar
    # Call below only for structured grid cases
    if not vmodel.attrs.get('mio_has_unstructured_grid',False):
        f = plt.gcf()
        model_ax = f.get_axes()[0]
        cax = f.get_axes()[1]
        #get the position of the plot axis and use this to rescale nicely the color bar to the height of the plot.
        position_m = model_ax.get_position()
        position_c = cax.get_position()
        cax.set_position([position_c.x0, position_m.y0, position_c.x1 - position_c.x0, (position_m.y1-position_m.y0)*1.1])
        cax.set_ylabel(ylabel,fontweight='bold',**text_kwargs)
        cax.tick_params(labelsize=text_kwargs['fontsize']*0.8,length=10.0,width=2.0,grid_linewidth=2.0)    
    
    #plt.tight_layout(pad=0)
    savefig(outname + '.png', loc=4, logo_height=100, dpi=150)
    return ax
    
def calculate_boxplot(df, column=None, label=None, plot_dict=None, comb_bx = None, label_bx = None):
    """Combines data into acceptable format for box-plot
    
    Parameters
    ----------
    df : dataframe
        Model/obs pair object
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
                 plot_dict=None, fig_dict=None,text_dict=None,debug=False):
    
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
    if debug == False:
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
        ylabel = label_bx[0]['column']
    
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
    savefig(outname + '.png', loc=4, logo_height=100)
