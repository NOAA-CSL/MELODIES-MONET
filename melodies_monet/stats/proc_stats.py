# Copyright (C) 2022 National Center for Atmospheric Research and National Oceanic and Atmospheric Administration
# SPDX-License-Identifier: Apache-2.0
#

#Simple MONET utility to calculate statistics from paired hdf file

import os
from glob import glob
import sys

import subprocess
from distutils.spawn import find_executable
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

import monet  
from monet.util.tools import calc_8hr_rolling_max,calc_24hr_ave,get_relhum
from monet.util.stats import STDO, STDP, MNB, MNE, MdnNB, MdnNE, NMdnGE, NO, NOP, NP, MO, MP, MdnO, MdnP, RM, RMdn, MB, MdnB, NMB, NMdnB, FB, ME, MdnE,NME, NMdnE, FE, MNPB, MdnNPB, MNPE, MdnNPE, NMPB, NMdnPB, NMPE, NMdnPE, R2, RMSE, d1, E1, IOA, AC, HSS, ETS, WDMB, WDMdnB, WDNMB_m, WDME, WDMdnE, WDRMSE, WDIOA, WDAC
import pandas as pd
import numpy as np
from numpy import sqrt
import matplotlib.pyplot as plt
from ..plots import savefig

def produce_stat_dict(stat_list,spaces=False):
    """Select statistics. Only statistics listed in the default dictionary 
    below are available.
    
    Parameters
    ----------
    stat_list : list of strings
        List of statistic abbreviations specified in input yaml file
    spaces : boolean
        Whether to leave spaces in the string containing the full name (True) 
        or remove spaces (False).
        
    Returns
    -------
    dictionary
        dictionary of statistics including abbreviations and full name

    """
    dict_stats_def = {'STDO' : 'Obs Standard Deviation', 
                      'STDP' : 'Mod Standard Deviation',
                      'MNB' : 'Mean Normalized Bias (%)',
                      'MNE' : 'Mean Normalized Gross Error (%)',
                      'MdnNB' : 'Median Normalized Bias (%)',
                      'MdnNE' : 'Median Normalized Gross Error (%)',
                      'NMdnGE' : 'Normalized Median Gross Error (%)',
                      'NO' : 'Obs Number',
                      'NOP' : 'Pairs Number',
                      'NP' : 'Mod Number',
                      'MO' : 'Obs Mean',
                      'MP' : 'Mod Mean',
                      'MdnO' : 'Obs Median',
                      'MdnP' : 'Mod Median',
                      'RM' : 'Mean Ratio Obs/Mod',
                      'RMdn' : 'Median Ratio Obs/Mod',
                      'MB' : 'Mean Bias',
                      'MdnB' : 'Median Bias',
                      'NMB' : 'Normalized Mean Bias (%)',
                      'NMdnB' : 'Normalized Median Bias (%)',
                      'FB' : 'Fractional Bias (%)',
                      'ME' : 'Mean Gross Error', 
                      'MdnE' : 'Median Gross Error', 
                      'NME' : 'Normalized Mean Error (%)',
                      'NMdnE' : 'Normalized Median Error (%)',
                      'FE' : 'Fractional Error (%)',
                      'R2' : 'Coefficient of Determination (R2)',
                      'RMSE' : 'Root Mean Square Error',
                      'd1' : 'Modified Index of Agreement',
                      'E1' : 'Modified Coefficient of Efficiency',
                      'IOA' : 'Index of Agreement',
                      'AC' : 'Anomaly Correlation'}
    stat_fullname_list = []
    for stat_id in stat_list:
        if spaces == False:
            stat_fullname_list.append(dict_stats_def[stat_id].replace(" ", "_"))
        else:
            stat_fullname_list.append(dict_stats_def[stat_id])
    #Note if you try to add a stat not in this list there will be an error in MELODIES-MONET, 
    #which is intended. If you want to add a new stat, please add the full name to the dictionary above.
    return stat_fullname_list

#Once redue the regulatory calculations. Use these from the surfplots. Or add these to a util script.
#Any stats not calculated in MONET will be added as a routine here.

def calc(df,stat=None,obsvar=None,modvar=None,wind=False):
    """Calculate statistics
    
    Parameters
    ----------
    df : dataframe
        model/obs pair data
    obsvar : str
        Column label of observation variable
    modvar : str
        Column label of model variable
    wind : boolean
        If variable is wind MONET applies a special calculation to handle 
        negative and positive values. If wind (True) if not wind (False).
        
    Returns
    -------
    real
        statistical value

    """
    obs = df[obsvar].values
    mod = df[modvar].values
   
    if stat == 'STDO':
        value = STDO(obs,mod,axis=None)
    elif stat == 'STDP':
        value = STDP(obs,mod,axis=None)
    #MNB looks wrong. Don't use for now. 
    elif stat == 'MNB':
        value = MNB(obs,mod,axis=None)
    #MNE looks wrong. Don't use for now. 
    elif stat == 'MNE':
        value = MNE(obs,mod,axis=None)
    elif stat == 'MdnNB':
        value = MdnNB(obs,mod,axis=None)
    elif stat == 'MdnNE':
        value = MdnNE(obs,mod,axis=None)
    elif stat == 'NMdnGE':
        value = NMdnGE(obs,mod,axis=None)
    elif stat == 'NO':
        value = NO(obs,mod,axis=None)
    elif stat == 'NOP':
        value = NOP(obs,mod,axis=None)
    elif stat == 'NP':
        value = NP(obs,mod,axis=None)
    elif stat == 'MO':
        value = MO(obs,mod,axis=None)
    elif stat == 'MP':
        value = MP(obs,mod,axis=None)
    elif stat == 'MdnO':
        value = MdnO(obs,mod,axis=None)
    elif stat == 'MdnP':
        value = MdnP(obs,mod,axis=None)
    elif stat == 'RM':
        value = RM(obs,mod,axis=None)
    elif stat == 'RMdn':
        value = RMdn(obs,mod,axis=None)
    elif stat == 'MB':
        if wind == True:
            value = WDMB(obs,mod,axis=None)
        else:
            value = MB(obs,mod,axis=None)
    elif stat == 'MdnB':
        if wind == True:
            value = WDMdnB(obs,mod,axis=None)
        else:
            value = MdnB(obs,mod,axis=None)
    elif stat == 'NMB':
        if wind == True:
            value = WDNMB_m(obs,mod,axis=None)
        else:
            value = NMB(obs,mod,axis=None)
    elif stat == 'NMdnB':
        value = NMdnB(obs,mod,axis=None)
    elif stat == 'FB':
        value = FB(obs,mod,axis=None)
    elif stat == 'ME':
        if wind == True:
            value = WDME(obs,mod,axis=None)
        else:
            value = ME(obs,mod,axis=None)
    elif stat == 'MdnE':
        if wind == True:
            value = WDMdnE(obs,mod,axis=None)
        else:
            value = MdnE(obs,mod,axis=None)
    elif stat == 'NME':
        value = NME(obs,mod,axis=None)
    elif stat == 'NMdnE':
        value = NMdnE(obs,mod,axis=None)
    elif stat == 'FE':
        value = FE(obs,mod,axis=None)
    elif stat == 'R2':
        value = R2(obs,mod,axis=None)
    elif stat == 'RMSE':
        if wind == True: 
            value = WDRMSE(obs,mod,axis=None)
        else:
            value = RMSE(obs,mod,axis=None)
    elif stat == 'd1':
        value = d1(obs,mod,axis=None)
    elif stat == 'E1':
        value = E1(obs,mod,axis=None)
    elif stat == 'IOA':
        if wind == True: 
            value = WDIOA(obs,mod,axis=None)
        else:
            value = IOA(obs,mod,axis=None)
    elif stat == 'AC':
        if wind == True:
            value = WDAC(obs,mod,axis=None)
        else:
            value = AC(obs,mod,axis=None)
    else:
        print('Stat not found: ' + stat)
        value = np.nan
            
    return value

def create_table(df,outname='plot',title='stats',out_table_kwargs=None,debug=False):
    """Calculates all of the specified statistics, save to csv file, and
    optionally save to a figure visualizing the table. 
    
    Parameters
    ----------
    df : dataframe
        model/obs pair data
    outname : str
        file location and name of plot (do not include .png)
    title : str
        Title to include on the figure visualizing the table
    out_table_kwargs : dictionary
        Dictionary containing information to create the figure visualizing the 
        table.
    debug : boolean
        Whether to plot interactively (True) or not (False). Flag for 
        submitting jobs to supercomputer turn off interactive mode.
        
    Returns
    -------
    csv file, plot
        csv file and optional plot containing the statistical calculations 
        specified in the input yaml file.

    """
    if debug == False:
        plt.ioff()
        
    #Define defaults if not provided:
    out_table_def = dict(fontsize=16.,xscale=1.2,yscale=1.2,figsize=[10,7],edges='open')
    if out_table_kwargs is not None:
        table_kwargs = {**out_table_def, **out_table_kwargs}
    else:
        table_kwargs = out_table_def
        
    #Create a table graphic
    fig, ax = plt.subplots(figsize=table_kwargs['figsize'])
    ax.axis('off')
    ax.axis('tight')
    
    rows=df['Stat_FullName'].values.tolist()
    
    df = df.drop(columns=['Stat_FullName'])
    
    t=ax.table(cellText=df.values, rowLabels=rows,
               colLabels=df.columns,loc='center',edges=table_kwargs['edges'])
    t.auto_set_font_size(False) 
    t.set_fontsize(table_kwargs['fontsize'])
    t.scale(table_kwargs['xscale'], table_kwargs['yscale'])
    plt.title(title,fontsize=table_kwargs['fontsize']*1.1,fontweight='bold')
    fig.tight_layout()
    savefig(outname + '.png', loc=1, logo_height=70)

    return


