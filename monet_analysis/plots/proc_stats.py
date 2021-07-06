#!/usr/bin/env python

###############################################################
# < next few lines under version control, D O  N O T  E D I T >
# $Date: 2018-03-29 10:12:00 -0400 (Thu, 29 Mar 2018) $
# $Revision: 100014 $
# $Author: Barry.Baker@noaa.gov $
# $Id: nemsio2nc4.py 100014 2018-03-29 14:12:00Z Barry.Baker@noaa.gov $
###############################################################

#Rebecca Schwantes adapted scripts from Patrick Campbell to MONET-analysis, added additional features, and generalized across all observations.

#Simple MONET utility to calculate statistics from paired hdf file

import os
from glob import glob
import sys

import subprocess
from distutils.spawn import find_executable
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

import monet  
from monet.util.tools import calc_8hr_rolling_max,calc_24hr_ave,get_relhum
from monet.util.stats import STDO, STDP, MNB, MNE, MdnNB, MdnNE, NMdnGE, NO, NOP, NP, MO, MP, MdnO, MdnP, RM, RMdn, MB, MdnB, NMB, NMdnB, FB, ME, MdnE,NME_m, NMdnE, FE, MNPB, MdnNPB, MNPE, MdnNPE, NMPB, NMdnPB, NMPE, NMdnPE, R2, RMSE, d1, E1, IOA_m, AC, HSS, ETS, WDMB_m, WDMdnB, WDNMB_m, WDME_m, WDMdnE, WDRMSE_m, WDIOA_m, WDAC
import pandas as pd
import numpy as np
from numpy import sqrt
import matplotlib.pyplot as plt
from new_models import code_to_move_to_monet as code_m_new

def produce_stat_dict(stat_list,spaces=False):
    #If spaces = True, leave spaces in string.
    #If spaces = False, remove spaces in string.
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
                      'MNPB' : 'Mean Normalized Peak Bias (%)',
                      'MdnNPB' : 'Median Normalized Peak Bias (%)',
                      'MNPE' : 'Mean Normalized Peak Error (%)',
                      'MdnNPE' : 'Median Normalized Peak Error (%)',
                      'NMPB' : 'Normalized Mean Peak Bias (%)',
                      'NMdnPB' : 'Normalized Median Peak Bias (%)',
                      'NMPE' : 'Normalized Mean Peak Error (%)',
                      'NMdnPE' : 'Normalized Median Peak Error (%)',
                      'R2' : 'Coefficient of Determination (R2)',
                      'RMSE' : 'Root Mean Square Error',
                      'd1' : 'Modified Index of Agreement',
                      'E1' : 'Modified Coefficient of Efficiency',
                      'IOA' : 'Index of Agreement',
                      'AC' : 'Anomaly Correlation',
                      'HSS' : 'Heidke Skill Score',
                      'ETS' : 'Equitable Threat Score'}
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

    obs = df[obsvar]
    mod = df[modvar]
   
    if stat == 'STDO':
        value = STDO(obs,mod,axis=None)
    elif stat == 'STDP':
        value = STDP(obs,mod,axis=None)
    #This looks wrong. Don't use for now. Will look into it.
    elif stat == 'MNB':
        value = MNB(obs,mod,axis=None)
    #MNE is erroring will look into soon. 
    #Need to figure out how this differs from NME?
    elif stat == 'MNE':
        value = MNE(obs,mod,axis=None)
    elif stat == 'MdnNB':
        value = MdnNB(obs,mod,axis=None)
    #MdnNE is also not working. Will look into this soon.
    elif stat == 'MdnNE':
        value = MdnNE(obs,mod,axis=None)
    #NMdnGE is also not working. Will look into this soon.
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
            value = WDMB_m(obs,mod,axis=None)
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
    #ME also not working, will look into this soon.
    elif stat == 'ME':
        if wind == True:
            value = WDME_m(obs,mod,axis=None)
        else:
            value = ME(obs,mod,axis=None)
    #MdnE also not working, will look into this soon.
    elif stat == 'MdnE':
        if wind == True:
            value = WDMdnE(obs,mod,axis=None)
        else:
            value = MdnE(obs,mod,axis=None)
    elif stat == 'NME':
        value = NME_m(obs,mod,axis=None)
    #NMdnE also not working, will look into this soon.
    #Likely need all of these to be converted to _m versions?
    elif stat == 'NMdnE':
        value = NMdnE(obs,mod,axis=None)
    #FE also not working.
    elif stat == 'FE':
        value = FE(obs,mod,axis=None)
    #For all peak stats, need to figure out what paxis is and add this
    #For now, these probably will not work.
    elif stat == 'MNPB':
        value = MNPB(obs,mod,axis=None)    
    elif stat == 'MdnNPB':
        value = MdnNPB(obs,mod,axis=None)
    elif stat == 'MNPE':
        value = MNPE(obs,mod,axis=None)
    elif stat == 'MdnNPE':
        value = MdnNPE(obs,mod,axis=None)
    elif stat == 'NMPB':
        value = NMPB(obs,mod,axis=None)
    elif stat == 'NMdnPB':
        value = NMdnPB(obs,mod,axis=None)
    elif stat == 'NMPE':
        value = NMPE(obs,mod,axis=None)
    elif stat == 'NMdnPE':
        value = NMdnPE(obs,mod,axis=None)
    #This ends the paxis types that I need to test more.
    elif stat == 'R2':
        value = R2(obs,mod,axis=None)
    elif stat == 'RMSE':
        if wind == True: 
            value = WDRMSE_m(obs,mod,axis=None)
        else:
            value = RMSE(obs,mod,axis=None)
    #Also does not work.
    elif stat == 'd1':
        value = d1(obs,mod,axis=None)
    #Also does not work.
    elif stat == 'E1':
        value = E1(obs,mod,axis=None)
    elif stat == 'IOA':
        if wind == True: 
            value = WDIOA_m(obs,mod,axis=None)
        else:
            value = IOA_m(obs,mod,axis=None)
    elif stat == 'AC':
        if wind == True:
            value = WDAC(obs,mod,axis=None)
        else:
            value = AC(obs,mod,axis=None)
    #In order to use HSS and ETS, need to specify min/max value.
    #Need to look into this one too.
    elif stat == 'HSS':
        value = HSS(obs,mod,axis=None)
    elif stat == 'ETS':
        value = ETS(obs,mod,axis=None)
    else:
        print('Stat not found')
        value = np.nan
            
    return value

def create_table(df,outname='plot',title='stats',out_table_kwargs=None):
    
    #Define defaults if not provided:
    out_table_def = dict(fontsize=16.,xscale=1.2,yscale=1.2,figsize=[10,6],edges='open')
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
    fig.tight_layout()
    plt.title(title,fontsize=table_kwargs['fontsize']*1.1,fontweight='bold')
    code_m_new.savefig(outname + '.png',loc=1, height=70, decorate=True, 
                                                   bbox_inches='tight', dpi=200)

    return


