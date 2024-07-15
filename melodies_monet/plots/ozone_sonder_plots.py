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
from matplotlib.colors import TwoSlopeNorm, ListedColormap, LinearSegmentedColormap
from matplotlib.patches import Rectangle
from matplotlib.ticker import FuncFormatter
from monet.util.tools import get_epa_region_bounds as get_epa_bounds
import math
from ..plots import savefig
from .surfplots import make_24hr_regulatory,calc_24hr_ave_v1,make_8hr_regulatory,calc_8hr_rolling_max_v1,calc_default_colors,new_color_map,map_projection,get_utcoffset,make_timeseries,make_taylor,calculate_boxplot,make_boxplot

#Define ozone sonder, vertical single date plot
def make_vertical_single_date(df,comb_bx,model_name_list,altitude_range,altitude_method,ozone_range,station_name,release_time,fill_color_list,plot_dict,fig_dict,text_dict):
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
    #set default figure size
    if fig_dict is not None:
        f,ax = plt.subplots(**fig_dict)
    else:
        f,ax = plt.subplots(figsize=(8,8))
    #set default text size
    def_text = dict(fontsize=20)
    if text_dict is not None:
        text_kwargs = {**def_text, **text_dict}
    else:
        text_kwargs = def_text

    #Make plot
    if altitude_method[0] =='sea level':
        ax.plot(O3_OBS,ALT_sl,'-*',color=fill_color_list[0],label = 'OBS:'+model_name_list[0])
        for i in range(len(O3_MODEL_ALL)): 
            ax.plot(O3_MODEL_ALL[i],ALT_sl,'-*',color=fill_color_list[i+1],label = 'MOD:'+model_name_list[i+1])
        plt.title('Comparison at '+str(station_name[0])+' on '+str(release_time)+' UTC')
        plt.legend()
        plt.ylim(altitude_range[0],altitude_range[1])
        plt.xlim(ozone_range[0],ozone_range[1])
        ax.set_ylabel('ALT-sea level (km)')
        ax.set_xlabel('O3 (ppbv)')
    elif altitude_method[0] == 'ground level':
        ax.plot(O3_OBS,ALT_gl,'-*',color=fill_color_list[0],label = 'OBS:'+model_name_list[0])
        for i in range(len(O3_MODEL_ALL)):
            ax.plot(O3_MODEL_ALL[i],ALT_gl, '-*',color=fill_color_list[i+1],label = 'MOD:'+model_name_list[i+1])     
        plt.title('Comparison at '+str(station_name[0])+' on '+str(release_time)+' UTC')
        plt.legend()
        plt.ylim(altitude_range[0],altitude_range[1])
        plt.xlim(ozone_range[0],ozone_range[1])
        plt.ylabel('ALT-ground level (km)')
        plt.xlabel('O3 (ppbv)')


#Define ozone sonder, vertical single date plot
def make_vertical_boxplot_os(df,comb_bx,label_bx,model_name_list,altitude_range,altitude_method,altitude_threshold_list,ozone_range,station_name,release_time,fill_color_list,plot_dict,fig_dict,text_dict):
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

    #set default figure size
    if fig_dict is not None:
        f,ax = plt.subplots(**fig_dict)
    else:
        f,ax = plt.subplots(figsize=(8,8))

    #set default text size
    def_text = dict(fontsize=20)
    if text_dict is not None:
        text_kwargs = {**def_text, **text_dict}
    else:
        text_kwargs = def_text

    #set color style
    obs_dict = dict(color='k', linestyle='-',marker='*', linewidth=2.0, markersize=5.)
    if plot_dict is not None:
        plot_kwargs = {**obs_dict, **plot_dict}
        plot_kwargs.pop('column')
    else:
        plot_kwargs = obs_dict

    #Define characteristics of boxplot.
    boxprops = {'edgecolor': 'k', 'facecolor': 'w','linewidth': 1.5}
    lineprops = {'color': 'k', 'linewidth': 1.5}
    boxplot_kwargs = {'boxprops': boxprops, 'medianprops': lineprops,
                  'whiskerprops':  lineprops, 'capprops': lineprops,
                  'widths':0.3, 'patch_artist': True,'showfliers': False
#                  'zorder':0,'showfliers': False
#                  'showmeans': False,'notch'= False,'showfliers': False
                  }
    
    #Make plot
    #get position list
    position_list_1 = altitude_threshold_list[:-1]
    position_list_2 = altitude_threshold_list[1:]
    position_list_mid = [(position_list_1[i]+position_list_2[i])/2 for i in range(len(position_list_1))]

    if altitude_method[0] =='sea level':
        output_list_obs = split_by_threshold(O3_OBS,ALT_sl,altitude_threshold_list)
        bplot_obs=ax.boxplot(output_list_obs,vert = False,**boxplot_kwargs,positions=position_list_mid)
        #add a line plot
        ax.plot(O3_OBS,ALT_sl,'-k*',label = 'OBS:'+model_name_list[0])

        for i in range(len(O3_MODEL_ALL)):
            output_list_model = split_by_threshold(O3_MODEL_ALL[i],ALT_sl,altitude_threshold_list)
            bplot_model=ax.boxplot(output_list_model,vert = False,**boxplot_kwargs,positions=position_list_mid)
            #add a line plot
            ax.plot(O3_MODEL_ALL[i],ALT_sl,'-*',color=fill_color_list[i+1],label = 'MODEL:'+model_name_list[i+1])

            colors=[fill_color_list[i+1]]*(len(altitude_threshold_list)-1)
            for patch,color in zip(bplot_model['boxes'],colors):
                patch.set_edgecolor(color)

        plt.title('Comparison at '+str(station_name[0])+' on '+str(release_time)+' UTC')
        plt.legend()
        plt.ylim(altitude_range[0],altitude_range[1])
        plt.xlim(ozone_range[0],ozone_range[1])
        plt.ylabel('ALT-sea level (km)')
        plt.xlabel('O3 (ppbv)')
    elif altitude_method[0] == 'ground level':
        output_list_obs = split_by_threshold(O3_OBS,ALT_gl,altitude_threshold_list)
        bplot_obs=ax.boxplot(output_list_obs,vert = False,**boxplot_kwargs,positions=position_list_mid)
        #add a line plot
        ax.plot(O3_OBS,ALT_gl,'-k*',label = 'OBS:'+model_name_list[0])
        for i in range(len(O3_MODEL_ALL)):
            output_list_model = split_by_threshold(O3_MODEL_ALL[i],ALT_gl,altitude_threshold_list)
            bplot_model=ax.boxplot(output_list_model,vert = False,**boxplot_kwargs,positions=position_list_mid)
            #add a line plot
            ax.plot(O3_MODEL_ALL[i],ALT_gl,'-*',color=fill_color_list[i+1],label = 'MODEL:'+model_name_list[i+1])
            
            colors=[fill_color_list[i+1]]*(len(altitude_threshold_list)-1)
            for patch,color in zip(bplot_model['boxes'],colors):
                patch.set_edgecolor(color)

        plt.title('Comparison at '+str(station_name[0])+' on '+str(release_time)+' UTC')
        plt.legend()
        plt.ylim(altitude_range[0],altitude_range[1])
        plt.xlim(ozone_range[0],ozone_range[1])
        plt.ylabel('ALT-ground level (km)')
        plt.xlabel('O3 (ppbv)')

def split_by_threshold(o3_list_input,alt_list_input,threshold_list_input):
    df = pd.DataFrame(data={'o3':o3_list_input,'alt':alt_list_input})
    output_list = []
    for i in range(1,len(threshold_list_input)):
        df_here = df.o3.loc[(df.alt>threshold_list_input[i-1])&(df.alt<=threshold_list_input[i])]
        output_list.append(df_here.values)
    return output_list

def density_scatter_plot_os(df,altitude_range,ozone_range,station_name,altitude_method,cmap_method):
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
    O3_OBS = df_short['o3']
    O3_MODEL = df_short['o3_ave']

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
    Modify_OBS.append(ozone_range[0])
    Modify_OBS.append(ozone_range[1])
    Modify_OBS.sort()
    plt.plot(Modify_OBS, slope(Modify_OBS),color='k',linestyle='-.',label='best fit')
   
    #plot Y=X line
    plt.axline((ozone_range[1],ozone_range[1]),slope=1,color='k',linestyle='-',label='Y=X')

    plt.xlim(ozone_range[0],ozone_range[1])
    plt.ylim(ozone_range[0],ozone_range[1])
    plt.xlabel('O$_3$ Obs (ppbv)')
    plt.ylabel('O$_3$ Model (ppbv)')
    plt.legend()




     








































