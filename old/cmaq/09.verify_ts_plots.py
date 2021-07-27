#!/usr/bin/env python

###############################################################
# < next few lines under version control, D O  N O T  E D I T >
# $Date: 2018-03-29 10:12:00 -0400 (Thu, 29 Mar 2018) $
# $Revision: 100014 $
# $Author: Barry.Baker@noaa.gov $
# $Id: nemsio2nc4.py 100014 2018-03-29 14:12:00Z Barry.Baker@noaa.gov $
###############################################################

__author__ = 'Patrick Campbell'
__email__ = 'Patrick.C.Campbell@noaa.gov'
__license__ = 'GPL'

import os
import subprocess
import sys
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser

import cartopy.crs as ccrs
import dask
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
import sys
import monet
from monet.util.tools import calc_8hr_rolling_max,calc_24hr_ave,get_relhum, kolmogorov_zurbenko_filter

sns.set_context('notebook')

plt.ioff()
'''
Simple utility to make time series plots
'''


def  make_24hr_regulatory(df,col=None):
     """ Make 24-hour averages """
     return calc_24hr_ave(df,col)

def  make_8hr_regulatory(df,col=None):
     """ Make 8-hour rolling average daily """
     return calc_8hr_rolling_max(df,col,window=8)

def make_kz_filter(df,col=None,window=None,iterations=None):
    """ Make kolmogorov_zurbenko_filter """
    return kolmogorov_zurbenko_filter(df,col,window,iterations)

def chdir(fname):
    dir_path = os.path.dirname(os.path.realpath(fname))
    os.chdir(dir_path)
    return os.path.basename(fname)


def load_paired_data(fname):
    return pd.read_hdf(fname)


def make_plots(df, variable, obs_variable, startdate, enddate, vmin, vmax, ylog, region, modcount, out_name):
    
    print(
            ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
    print('Creating Plot:', obs_variable, 'for period:', startdate, 'to ', enddate  )
    print(
            ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")    
 
    make_timeseries_epa(df, out_name, startdate, enddate, vmin, vmax, ylog, region, modcount,col1=obs_variable, col2=variable)


def make_timeseries_epa(
        df,
        savename,
        startdate,
        enddate,
        vmin,
        vmax,
        ylog,
        region,
        modcount,
        col1='OZONE',
        col2='O3'
):
    from monet.util.tools import get_epa_region_df as epard
    from monet.plots import savefig
    import seaborn as sns

    df.index=df.time
    print(subset_name)
    if subset_name == 'epa_region':
     df.query('epa_region == '+'"'+region+'"',inplace=True)
    if subset_name == 'state_name':
     df.query('state_name == '+'"'+region+'"',inplace=True)
    if subset_name == 'siteid':
     print('subsetting...')
     df.query('siteid == '+'"'+region+'"',inplace=True)
    if modcount == 0:
     ax=df[col1].resample('H').mean().plot(marker='.',color='darkslategrey',label='OBS')
     ax=df[col2].resample('H').mean().plot(ax=ax,label='MOD')   
    if modcount == 1:
     ax=df[col2].resample('H').mean().plot(color='red',label='MOD2')
    if modcount == 2:
     ax=df[col2].resample('H').mean().plot(color='green',label='MOD3')
    if modcount == 3:
     ax=df[col2].resample('H').mean().plot(color='purple',label='MOD4')
    if modcount == 4:
     ax=df[col2].resample('H').mean().plot(color='orange',label='MOD5')
    sns.despine()
    if vmin != None and vmax != None:
     plt.ylim(vmin, vmax)
    if ylog == True:
     plt.yscale('log')
    plt.xlim([startdate,enddate])
    plt.xlabel('TIME (UTC)')
    plt.ylabel(r'Concentration (ppbV or micrograms/m-3)')
    plt.legend(loc=1)
    plt.tight_layout(pad=0)
    name = "{}.ts.pdf".format(savename)
    monet.plots.savefig(name, dpi=100, loc=4, decorate=False)

if __name__ == '__main__':

    parser = ArgumentParser(
        description='Make Time Series Plots for each time step in files',
        formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '-p',
        '--paired_data',
        help='paired data input file names',
        type=str,
        nargs='+',
        required=True)
    parser.add_argument(
        '-s', '--species', nargs='+', help='Species', required=False, default=['OZONE'])
    parser.add_argument(
        '-n',
        '--output_name',
        help='Box-whisker plot Output base name',
        type=str,
        required=False,
        default='CMAQ_AIRNOW')
    parser.add_argument(
        '-r',
        '--regulatory',
        help='boolean set to True fore 8-hrmax  or 24-ave NAAQS regulatory calcs',
        type=bool,
        required=False,
        default=False)
    parser.add_argument(
        '-kz',
        '--kzfilter',
        help='boolean set to True for applying a kolmogorov_zurbenko_filter',
        type=bool,
        required=False,
        default=False)
    parser.add_argument(
        '-sd',
        '--startdate',
        help='Startdate for time series over a period YYYY-MM-DD HH:MM:SS',
        type=str,
        required=False,
        default=None)
    parser.add_argument(
        '-ed',
        '--enddate',
        help='Enddate for time series over a period YYYY-MM-DD HH:MM:SS',
        type=str,
        required=False,
        default=None)
    parser.add_argument(
        '-e',
        '--epa_region',
        help='EPA Region acronymn, state name, or siteid',
        required=False,
        default='R1')
    parser.add_argument(
        '-sn',
        '--subset_name',
        help='name of subset type (epa_region, state_name, or siteid)',
        required=False,
        default='epa_region')
    parser.add_argument(
        '-miny', '--miny_scale', help='Set static min y-scale', type=float, required=False, default=None)
    parser.add_argument(
        '-maxy', '--maxy_scale', help='Set static max y-scale', type=float, required=False, default=None)
    parser.add_argument(
        '-ylog', '--ylog_scale', help='Set log y-scale', type=bool, required=False, default=False)
    parser.add_argument(
        '-kzw', '--kz_window', help='Set window, m, in model time units  (default = 15 days (360 hrs), AQ Baseline Component, Weiss and Comrie, (2005))', 
         type=int, required=False, default=360)
    parser.add_argument(
        '-kzi', '--kz_iterations', help='Set kz iterations (default = 5, AQ Baseline Component, Weiss and Comrie (2005), removes cycles < 33 days)', 
         type=int, required=False, default=5)

    args = parser.parse_args()

    paired_data = args.paired_data
    species     = args.species
    out_name    = args.output_name
    startdate   = args.startdate
    enddate     = args.enddate
    reg         = args.regulatory
    region      = args.epa_region
    subset_name = args.subset_name
    vmin        = args.miny_scale
    vmax        = args.maxy_scale
    ylog        = args.ylog_scale
    kzfilt      = args.kzfilter
    kzw         = args.kz_window
    kzi         = args.kz_iterations

    for jj in species:
     plt.close() 
     modcount = 0
     for xx in paired_data:
      print('model#')
      print((modcount+1))  
#load the paired dataframe
      df = load_paired_data(xx)
      mapping_table = {'OZONE':'O3', 'PM2.5':'PM25_TOT', 'PM10':'PM10_new', 'CO':'CO_new', 'NO':'NO_new', 'NO2':'NO2_new', 'SO2':'SO2_new','NOX':'NOX_new','NOY':'NOY_new','TEMP':'TEMP2','WS':'WSPD10','WD':'WDIR10','SRAD':'GSW','BARPR':'PRSFC','PRECIP':'RT','RHUM':'Q2'}
      sub_map = {i: mapping_table[i] for i in species if i in mapping_table}
  
      df_replace = df.replace(0.0,np.nan) #Replace all exact 0.0 values with nan
      df_drop=df_replace.dropna(subset=[jj,sub_map.get(jj)]) #Drops all rows with obs species = NaN

#Converts OZONE, PM10, or PM2.5 dataframe to NAAQS regulatory values
      if jj == 'OZONE' and reg is True:
       df2 = make_8hr_regulatory(df_drop,[jj,sub_map.get(jj)]).rename(index=str,columns={jj+'_y':jj,sub_map.get(jj)+'_y':sub_map.get(jj)})
      elif jj == 'PM2.5' and reg is True:
       df2 = make_24hr_regulatory(df_drop,[jj,sub_map.get(jj)]).rename(index=str,columns={jj+'_y':jj,sub_map.get(jj)+'_y':sub_map.get(jj)})
      elif jj == 'PM10' and reg is True:
       df2 = make_24hr_regulatory(df_drop,[jj,sub_map.get(jj)]).rename(index=str,columns={jj+'_y':jj,sub_map.get(jj)+'_y':sub_map.get(jj)})
      else:
       df2=df_drop
#Convert airnow met variable if necessary:
      if jj == 'WS':
       df2.loc[:,'WS']=df2.loc[:,'WS']*0.514  #convert obs knots-->m/s
       df2.query('WS > 0.2',inplace=True)  #Filter out calm WS obs (< 0.2 m/s), should not be trusted--creates artificially larger postive  model bias
      elif jj == 'BARPR':
       df2.loc[:,'PRSFC']=df2.loc[:,'PRSFC']*0.01 #convert model Pascals-->millibars
      elif jj == 'PRECIP':
       df2.loc[:,'PRECIP']=df2.loc[:,'PRECIP']*0.1 #convert obs mm-->cm
      elif jj == 'TEMP':
       df2.loc[:,'TEMP2'] = df2.loc[:,'TEMP2']-273.16 #convert model K-->C
      elif jj == 'RHUM':
      #convert model mixing ratio to relative humidity
       df2.loc[:,'Q2'] = get_relhum(df2.loc[:,'TEMP2'],df2.loc[:,'PRSFC'],df2.loc[:,'Q2'])
      elif jj == 'CO':
       df2.loc[:,'CO']=df2.loc[:,'CO']*1000.0 #convert obs ppm-->ppb
      else:
       df2=df2
#subset for period, or use output frequency
      if startdate != None and enddate != None:

       mask = (df2['time'] >= startdate) & (df2['time'] <= enddate)
       dfnew =df2.loc[mask]
       import datetime
       startdatename_obj = datetime.datetime.strptime(startdate, '%Y-%m-%d %H:%M:%S')
       enddatename_obj   = datetime.datetime.strptime(enddate, '%Y-%m-%d %H:%M:%S')
       startdatename = str(datetime.datetime.strftime(startdatename_obj,'%Y-%m-%d_%H'))
       enddatename = str(datetime.datetime.strftime(enddatename_obj,'%Y-%m-%d_%H'))
       outname = "{}.{}.{}.{}.{}".format(out_name, region, jj, startdatename, enddatename)
       if reg is True:
        outname = "{}.{}.{}.{}.{}.{}".format(out_name, region, jj,startdatename, enddatename,'reg')
       if kzfilt is True:
        outname = "{}.{}.{}.{}.{}.{}".format(out_name, region, jj,startdatename, enddatename,'kz')
       if kzfilt is True and reg is True:
        outname = "{}.{}.{}.{}.{}.{}.{}".format(out_name, region, jj,startdatename, enddatename,'kz', 'reg')
       if jj == 'PM2.5':
        outname = outname.replace('PM2.5','PM2P5')
       if region == 'domain':
        outname = outname.replace('domain','5X')
      else:
       dfnew = df2
       outname = "{}.{}.{}".format(out_name,region, jj)
       if reg is True:
        outname = "{}.{}.{}.{}".format(out_name,region, jj, 'reg')
       if kzfilt is True:
        outname = "{}.{}.{}.{}".format(out_name,region, jj, 'kz')
       if kzfilt is True and reg is True:
        outname = "{}.{}.{}.{}.{}".format(out_name,region, jj, 'kz', 'reg')
       if jj == 'PM2.5':
        outname = outname.replace('PM2.5','PM2P5')
       if region == 'domain':
        outname = outname.replace('domain','5X')

      dfnew_drop=dfnew.dropna(subset=[jj,sub_map.get(jj)])

      initial_datetime = dfnew_drop.time.min()
#Appling a low-pass K-Z filter
      if kzfilt is True:
       dfnew_drop2 = make_kz_filter(dfnew_drop,[jj,sub_map.get(jj)],window=kzw,iterations=kzi).rename(index=str,columns={jj+'_y':jj,sub_map.get(jj)+'_y':sub_map.get(jj)})
      else:
       dfnew_drop2 = dfnew_drop
# make the plots
      make_plots(dfnew_drop2, sub_map.get(jj), jj, startdate, enddate, vmin, vmax, ylog, region, modcount,outname)
      modcount=modcount+1
