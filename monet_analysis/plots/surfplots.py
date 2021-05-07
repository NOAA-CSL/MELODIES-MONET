#!/usr/bin/env python

###############################################################
# < next few lines under version control, D O  N O T  E D I T >
# $Date: 2018-03-29 10:12:00 -0400 (Thu, 29 Mar 2018) $
# $Revision: 100014 $
# $Author: Barry.Baker@noaa.gov $
# $Id: nemsio2nc4.py 100014 2018-03-29 14:12:00Z Barry.Baker@noaa.gov $
###############################################################

#Original scripts by Patrick Campbell. Adapted to MONET-analysis by Rebecca Schwantes and Barry Baker

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

# from util import write_ncf

def make_24hr_regulatory(df, col=None):
    """ Make 24-hour averages """
    return calc_24hr_ave(df, col)


def make_8hr_regulatory(df, col=None):
    """ Make 8-hour rolling average daily """
    return calc_8hr_rolling_max(df, col, window=8)

def calc_default_colors(p_index):
    """ Use default colors """
    x = [dict(color='b', linestyle='--',marker='x'),
         dict(color='g', linestyle='-.',marker='o'),
         dict(color='r', linestyle=':',marker='v'),
         dict(color='c', linestyle='--',marker='^'),
         dict(color='m', linestyle='-.',marker='s')]
    #Repeat these 5 instances over and over if more than 5 lines.
    return x[p_index % 5]

def new_color_map():
    top = mpl.cm.get_cmap('Blues_r', 128)
    bottom = mpl.cm.get_cmap('Oranges', 128)
    newcolors = np.vstack((top(np.linspace(0, 1, 128)),
                           bottom(np.linspace(0, 1, 128))))
    return ListedColormap(newcolors, name='OrangeBlue')

def get_df_region(obj, region):
    from monet.util.tools import get_epa_region_df as get_epa
    if region.lower() == 'domain':
        obj['EPA_ACRO'] = 'domain'
        return obj
    else:
        obj = get_epa(region)
        return obj.loc[obj.EPA_ACRO == region.upper()]

def make_spatial_bias(df, column_o=None, label_o=None, column_m=None, 
                      label_m=None, ylabel = None, outname = 'plot', vmin=None, vmax=None, 
                      fig_dict=None, text_dict=None, map_dict=None):
    """Creates the MONET-Analysis spatial bias plot."""
    #Need to add possibility of including map_kwargs into yaml file for now use default.
    #Also need to add possibility of val_min and val_max this has to be in the obs though because dependent on species.
    def_map = dict(states=True)
    if map_dict is not None:
        map_kwargs = {**def_map, **map_dict}
    else:
        map_kwargs = def_map
  
    #set default text size
    def_text = dict(fontsize=14)
    if text_dict is not None:
        text_kwargs = {**def_text, **text_dict}
    else:
        text_kwargs = def_text
    
    cmap = new_color_map()
    #Take the mean for each siteid
    df_mean=df.groupby(['siteid'],as_index=False).mean()
    
    ax = monet.plots.sp_scatter_bias(
        df, col1=column_o, col2=column_m, map_kwargs=map_kwargs,val_max=vmax,val_min=vmin,
                    cmap=cmap, edgecolor='k',linewidth=.8)
    #date = pd.Timestamp(date)
    #dt = date - initial_datetime
    #dtstr = str(dt.days * 24 + dt.seconds // 3600).zfill(3)
    #plt.title(date.strftime('time=%Y/%m/%d %H:00 | CMAQ - AIRNOW '))
        
    #if region == 'domain':
    latmin= 25.0
    lonmin=-130.0
    latmax= 55.0
    lonmax=-55.0
    #else:
    # from monet.util.tools import get_epa_region_bounds as get_epa_bounds    
    # latmin,lonmin,latmax,lonmax,acro = get_epa_bounds(index=None,acronym=region)
   
    plt.xlim([lonmin,lonmax])
    plt.ylim([latmin,latmax]) 
  
    plt.tight_layout(pad=0)
    monet.plots.savefig(outname + '.png', bbox_inches='tight', dpi=200, decorate=True)
    
def make_timeseries(df, column=None, label=None, ax=None, avg_window=None, ylabel=None,
                    plot_dict=None, fig_dict=None, text_dict=None):
    """Creates the MONET-Analysis time series plot.

    Parameters
    ----------
    df : type
        Description of parameter `df`.
    label : type
        Description of parameter `label`.
    plot_dict : type
        Description of parameter `plot_dict`.

    Returns
    -------
    type
        Description of returned object.

    """
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
    
    #Then, if no plot has been created yet, create a plot and plot the obs.
    if ax is None: 
        #First define the colors for the observations.
        obs_dict = dict(color='k', linestyle='-',marker='x', linewidth=1.2, markersize=6.)
        if plot_dict is not None:
            #Whatever is not defined in the yaml file is filled in with the obs_dict here.
            plot_kwargs = {**obs_dict, **plot_dict}
        else:
            plot_kwargs = obs_dict
        if label is not None:
            plot_kwargs['label'] = label
        #scale the fontsize for the x and y labels by the text_kwargs
        plot_kwargs['fontsize'] = text_dict['fontsize']*0.8
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
        if label is not None:
            plot_dict['label'] = label
        #scale the fontsize for the x and y labels by the text_kwargs
        plot_dict['fontsize'] = text_dict['fontsize']*0.8
        if avg_window is None:
            ax = df[column].plot(ax=ax, **plot_dict)
        else:
            ax = df[column].resample(avg_window).mean().plot(ax=ax, legend=True, **plot_dict)    
    
    #Set parameters for all plots
    ax.set_ylabel(ylabel,fontweight='bold',**text_kwargs)
    ax.set_xlabel(df.index.name,fontweight='bold',**text_kwargs)
    ax.legend(frameon=False,fontsize=text_kwargs['fontsize']*0.8)
    ax.tick_params(axis='both',length=10.0,direction='inout')
    ax.tick_params(axis='both',which='minor',length=5.0,direction='out')
    return ax
    
def make_taylor(df, column_o=None, label_o='Obs', column_m=None, label_m='Model', dia=None, ylabel=None, ty_scale=1.5,
                    plot_dict=None, fig_dict=None, text_dict=None):
    """Creates the MONET-Analysis taylor plot.

    Parameters
    ----------
    df : type
        Description of parameter `df`.
    label : type
        Description of parameter `label`.
    plot_dict : type
        Description of parameter `plot_dict`.

    Returns
    -------
    type
        Description of returned object.

    """
    #First define items for all plots
    
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
    plt.legend(frameon=False,fontsize=text_kwargs['fontsize']*0.8)
    plt.xlabel('Standard deviation: '+ylabel,**text_kwargs)
    plt.ylabel('Standard deviation: '+ylabel,**text_kwargs)
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
    #plt.xticks(**text_kwargs)
    #plt.yticks(**text_kwargs)
    return dia

def make_spatial_overlay(paired, plot_dict=None):
    # TODO: write wrapper for overlay plots
        a = 1
