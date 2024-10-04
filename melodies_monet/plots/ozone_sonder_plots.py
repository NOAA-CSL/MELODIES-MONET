# Copyright (C) 2022 National Center for Atmospheric Research and National Oceanic and Atmospheric Administration
# SPDX-License-Identifier: Apache-2.0
#

import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
sns.set_context('paper')

#Define ozone sonder, vertical single date plot
def make_vertical_single_date(df,comb_bx,altitude_range,altitude_method,vmin, vmax,station_name,release_time,label_bx,fig_dict,text_dict):
    ALT_sl = df['altitude']
    O3_OBS = comb_bx[comb_bx.columns[0]].to_list()
    len_combx = np.shape(comb_bx)[1]
    O3_MODEL_ALL = []
    for i in range(1,len_combx):
        O3_MODEL = comb_bx[comb_bx.columns[i]].to_list()
        O3_MODEL_ALL.append(O3_MODEL)
    #release height info,get height of each release site to be substract    
    df_height = pd.DataFrame({
         'station':['Boulder, Colorado','Huntsville, Alabama','University of Rhode Island','Trinidad Head, California'],
         'height':[1.743,0.203,0.021,0.046]
         })
    if altitude_method[0] == 'ground level':
        height_value = df_height.loc[df_height['station']==station_name[0]].values[0][1]
        ALT_gl = ALT_sl - height_value
    #1 set default figure size
    if fig_dict is not None:
        f,ax = plt.subplots(**fig_dict)
    else:
        f,ax = plt.subplots(figsize=(8,8))
    #2 set default text size
    def_text = dict(fontsize=20)
    if text_dict is not None:
        text_kwargs = {**def_text, **text_dict}
    else:
        text_kwargs = def_text
    #3 set default plot key words
    ylabel = label_bx[0]['column']
    for i in range(len(label_bx)):
        label_bx[i].pop('column')

#Make plot
    if altitude_method[0] =='sea level':
        alt_p = ALT_sl
        alt_p_name = 'ALT-sea level (km)'
    elif altitude_method[0] == 'ground level':
        alt_p = ALT_gl
        alt_p_name = 'ALT-ground level (km)'

    plot_dict = label_bx[0]
    ax.plot(O3_OBS,alt_p,**plot_dict)
    for i in range(len(O3_MODEL_ALL)):
        plot_dict = label_bx[i+1]
        ax.plot(O3_MODEL_ALL[i],alt_p,**plot_dict)
    plt.title('Comparison at '+str(station_name[0])+' on '+str(release_time)+' UTC',**text_kwargs)
    plt.legend()
    plt.ylim(altitude_range[0],altitude_range[1])
    plt.xlim(vmin,vmax)
    ax.set_ylabel(alt_p_name,fontsize=text_kwargs['fontsize']*0.8)
    ax.set_xlabel(ylabel+' (ppbv)',fontsize=text_kwargs['fontsize']*0.8)

#Define ozone sonder, vertical single date plot
def make_vertical_boxplot_os(df,comb_bx,label_bx,altitude_range,altitude_method,vmin, vmax,altitude_threshold_list,station_name,release_time,fig_dict,text_dict):
    ALT_sl = df['altitude']
    O3_OBS = comb_bx[comb_bx.columns[0]].to_list()
    len_combx = np.shape(comb_bx)[1]
    O3_MODEL_ALL = []

    for i in range(1,len_combx):
        O3_MODEL = comb_bx[comb_bx.columns[i]].to_list()
        O3_MODEL_ALL.append(O3_MODEL)

    #release height info,get height of each release site to be substract 
    df_height = pd.DataFrame({
         'station':['Boulder, Colorado','Huntsville, Alabama','University of Rhode Island','Trinidad Head, California'],
         'height':[1.743,0.203,0.021,0.046]
         })
    if altitude_method[0] == 'ground level':
        height_value = df_height.loc[df_height['station']==station_name[0]].values[0][1]
        ALT_gl = ALT_sl - height_value

    #1 set default figure size
    if fig_dict is not None:
        f,ax = plt.subplots(**fig_dict)
    else:
        f,ax = plt.subplots(figsize=(8,8))

    #2 set default text size
    def_text = dict(fontsize=20)
    if text_dict is not None:
        text_kwargs = {**def_text, **text_dict}
    else:
        text_kwargs = def_text

    #3 set color style
    ylabel = label_bx[0]['column']
    for i in range(len(label_bx)):
        label_bx[i].pop('column')
    
    #Make plot
    #get position list
    position_list_1 = altitude_threshold_list[:-1]
    position_list_2 = altitude_threshold_list[1:]
    position_list_mid = [(position_list_1[i]+position_list_2[i])/2 for i in range(len(position_list_1))]

    if altitude_method[0] =='sea level':
        alt_p = ALT_sl
        alt_p_name = 'ALT-sea level (km)'
    elif altitude_method[0] == 'ground level':
        alt_p = ALT_gl
        alt_p_name = 'ALT-ground level (km)'

    output_list_obs = split_by_threshold(O3_OBS,alt_p,altitude_threshold_list)
    bplot_obs=ax.boxplot(output_list_obs,vert = False,patch_artist=True,
                         whiskerprops=dict(color=label_bx[0]['color']),
                         capprops=dict(color=label_bx[0]['color']),
                         boxprops=dict(facecolor='w',color=label_bx[0]['color']),
                         medianprops=dict(color=label_bx[0]['color']),
                         positions=position_list_mid)
    #add a line plot
    plot_dict = label_bx[0]
    ax.plot(O3_OBS,alt_p,**plot_dict)

    for i in range(len(O3_MODEL_ALL)):
        output_list_model = split_by_threshold(O3_MODEL_ALL[i],alt_p,altitude_threshold_list)
        bplot_model=ax.boxplot(output_list_model,vert = False,patch_artist=True,
                               whiskerprops=dict(color=label_bx[i+1]['color']),
                               capprops=dict(color=label_bx[i+1]['color']),
                               boxprops=dict(facecolor='w',color=label_bx[i+1]['color']),
                               medianprops=dict(color=label_bx[i+1]['color']),
                               positions=position_list_mid)
        #add a line plot
        plot_dict = label_bx[i+1]
        ax.plot(O3_MODEL_ALL[i],alt_p,**plot_dict)

    plt.title('Comparison at '+str(station_name[0])+' on '+str(release_time)+' UTC',**text_kwargs)
    plt.legend()
    plt.ylim(altitude_range[0],altitude_range[1])
    plt.xlim(vmin,vmax)
    plt.ylabel(alt_p_name,fontsize=text_kwargs['fontsize']*0.8)
    plt.xlabel(ylabel+' (ppbv)',fontsize=text_kwargs['fontsize']*0.8)

def split_by_threshold(o3_list_input,alt_list_input,threshold_list_input):
    df = pd.DataFrame(data={'o3':o3_list_input,'alt':alt_list_input})
    output_list = []
    for i in range(1,len(threshold_list_input)):
        df_here = df.o3.loc[(df.alt>threshold_list_input[i-1])&(df.alt<=threshold_list_input[i])]
        output_list.append(df_here.values)
    return output_list


def density_scatter_plot_os(df,altitude_range,vmin,vmax,station_name,altitude_method,cmap_method,modvar,obsvar):
    #release height info,get height of each release site to be substract    
    df_height = pd.DataFrame({
         'station':['Boulder, Colorado','Huntsville, Alabama','University of Rhode Island','Trinidad Head, California'],
         'height':[1.743,0.203,0.021,0.046]
         })
    if altitude_method[0] == 'ground level':
        height_value = df_height.loc[df_height['station']==station_name[0]].values[0][1]
    elif altitude_method[0] == 'sea level':
        height_value = 0

    #get o3 model, o3 sonder (obs) and height
    df_short = df[df['altitude']<altitude_range[1]+height_value]
    ALT = df_short['altitude']-height_value
    O3_OBS = df_short[obsvar]  #'o3'
    O3_MODEL = df_short[modvar]  #'o3_ave'

    #plot scatter and colorbar
    sc=plt.scatter(O3_OBS,O3_MODEL,c= ALT,vmin=altitude_range[0],vmax=altitude_range[1],cmap = 'jet',edgecolors='k',linewidth=0.5,s = 15)
    cb=plt.colorbar(sc)
    if altitude_method[0] == 'ground level':
        cb.set_label('Ground Level Altitude (km)',fontsize=15)
    elif altitude_method[0] == 'sea level':
        cb.set_label('Sea Level Altitude (km)',fontsize=15)

    #add some points to make best fit go entire domain
    slope = np.poly1d(np.polyfit(O3_OBS, O3_MODEL, 1))
    Modify_OBS = O3_OBS.to_list()
    Modify_OBS.append(vmin)
    Modify_OBS.append(vmax)
    Modify_OBS.sort()
    plt.plot(Modify_OBS, slope(Modify_OBS),color='k',linestyle='-.',label='best fit')
   
    #plot Y=X line
    plt.axline((vmax,vmax),slope=1,color='k',linestyle='-',label='Y=X')
    plt.xlim(vmin,vmax)
    plt.ylim(vmin,vmax)
    plt.xlabel(obsvar+' Obs (ppbv)')
    plt.ylabel(obsvar+' Model (ppbv)')
    plt.legend()
