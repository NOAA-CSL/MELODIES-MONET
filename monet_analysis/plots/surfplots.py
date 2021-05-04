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
sns.set_context('paper')

# from util import write_ncf

def make_24hr_regulatory(df, col=None):
    """ Make 24-hour averages """
    return calc_24hr_ave(df, col)


def make_8hr_regulatory(df, col=None):
    """ Make 8-hour rolling average daily """
    return calc_8hr_rolling_max(df, col, window=8)

def make_taylor_plot(df, model_var, obs_var , outname, obs_label, model_label, plot_dict=None, region=None,
                     epa_regulatory=False,time_avg=False):
    #If region is true query for only that specific region
    if region != None:
        df.query('epa_region == ' + '"' + region + '"', inplace=True)
    if time_avg == False:
        for t in df.time.unique():
            date = pd.Timestamp(t)
            print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
            print('Creating Plot:', obs_var, 'at time:', date)
            print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
            odf = df.loc[df.time == date, ['time', 'latitude', 'longitude', obs_var, model_var]]
                #add the time at the end of the file
            out_name = "{}.{}".format(outname, date)
            if ~odf.empty:
                make_taylor_diagram(odf, col1=obs_var, col2=model_var, label1=obs_label, label2=model_label,
                                        savename=out_name, plot_dict=plot_dict)
    # make total period taylor plot
    else:
        print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
        print('Creating Plot:', obs_var, 'for whole time period')
        print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
        make_taylor_diagram(df, col1=obs_var, col2=model_var, label1=obs_label, label2=model_label,
                            savename=outname, plot_dict=plot_dict)

def make_taylor_diagram(df, col1, col2, label1, label2, savename, plot_dict=None):
    #If we want to add fig_height and fig_width we will need to change plots.py in MONET itself.
    if plot_dict != None: #If specifiy ploting commands use them
        dia = monet.plots.plots.taylordiagram(df, col1=col1, col2=col2, label1=label1, label2=label2,
                                                  scale=plot_dict['taylor_diagram_scale'])
    else: #Else use defaults
        dia = monet.plots.plots.taylordiagram(df, col1=col1, col2=col2, label1=label1, label2=label2)
    plt.legend(loc=(0.8, 0.8))
    name = "{}.png".format(savename)
    monet.plots.savefig(name, bbox_inches='tight', dpi=100, loc=3, decorate=False)
    return dia
#Because this returns dia if we want you could update this to plot more than one model on a graph.

def make_spatial_bias(paired, plot_dict=None):
    # TODO: create wrapper for spatial bias
    a = 1

def make_timeseries(df, column=None, label=None, ax=None, plot_dict=None, fig_dict=None):
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
    # assume that 'time' is in the DataFrame
    if plot_dict is None: # assume that no plot has been created
        default_dict = dict(color='k',linewidth=1.2, linestyle='-',marker='x')
        if label is not None:
            default_dict['label'] = label
        # create the figure
        if fig_dict is not None:
            f,ax = plt.subplots(**fig_dict)    
        else: 
            f,ax = plt.subplots(figsize=(10,6))
        # plot the line
        ax = df.plot(ax=ax,x='time_local',y=column, **default_dict)
        ax.set_ylabel = column
    elif 'ax' != None:
        # this means that an axis handle already exists and use it to plot onto
        if label is not None:
            plot_dict['label'] = label
        ax = df.plot(ax=ax, x='time_local',y=column, **plot_dict)
    return ax
    

def make_spatial_overlay(paired, plot_dict=None):
    # TODO: write wrapper for overlay plots
        a = 1
