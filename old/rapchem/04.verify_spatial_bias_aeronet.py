#!/usr/bin/env python

###############################################################
# < next few lines under version control, D O  N O T  E D I T >
# $Date: 2018-03-29 10:12:00 -0400 (Thu, 29 Mar 2018) $
# $Revision: 100014 $
# $Author: Barry.Baker@noaa.gov $
# $Id: nemsio2nc4.py 100014 2018-03-29 14:12:00Z Barry.Baker@noaa.gov $
###############################################################

__author__ = 'Patrick.C.Campbell'
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
import seaborn as sns
import numpy as np
import monet
from monet.util.tools import calc_8hr_rolling_max,calc_24hr_ave,get_relhum
import dask.dataframe as dd

sns.set_context('notebook')

plt.ioff()
'''
Simple utility to make spatial plots from the NAQFC forecast and overlay observations
'''

initial_datetime = None

def  make_24hr_regulatory(df,col=None):
     """ Make 24-hour averages """
     return calc_24hr_ave(df,col)

def  make_8hr_regulatory(df,col=None):
     """ Make 8-hour rolling average daily """
     return calc_8hr_rolling_max(df,col,window=8)


def chdir(fname):
    dir_path = os.path.dirname(os.path.realpath(fname))
    os.chdir(dir_path)
    return os.path.basename(fname)

def load_paired_data(fname):
    return dd.read_hdf(fname,'/*').compute()


def make_spatial_bias_plot(df,
                           out_name,
                           vmin,
                           vmax,
                           col1='aod_550nm',
                           col2='pm25aod550',
                           date=None,
                           region='domain',
                           **kwargs):
    if region == 'domain':
     ax = monet.plots.sp_scatter_bias(
        df, col1=col1, col2=col2, map_kwargs=dict(states=False),val_max=vmax,val_min=vmin,**kwargs)
    else:
     ax = monet.plots.sp_scatter_bias(
        df, col1=col1, col2=col2, map_kwargs=dict(states=True),val_max=vmax,val_min=vmin,**kwargs)

    date = pd.Timestamp(date)
    dt = date - initial_datetime
    dtstr = str(dt.days * 24 + dt.seconds // 3600).zfill(3)
    plt.title(date.strftime('time=%Y/%m/%d %H:00 | CMAQ - AIRNOW '))
        
    if region == 'domain':
     latmin=-90.0
     lonmin=-180.0
     latmax=90.0
     lonmax=180.0
    else:
     from monet.util.tools import get_giorgi_region_bounds as get_giorgi_bounds    
     latmin,lonmin,latmax,lonmax,acro = get_giorgi_bounds(index=None,acronym=region)
   
    plt.xlim([lonmin,lonmax])
    plt.ylim([latmin,latmax]) 
  
    plt.tight_layout(pad=0)
    savename = "{}.{}.{}.jpg".format(out_name,
                                     initial_datetime.strftime('spbias'),
                                     dtstr)
    print(savename)
    monet.plots.savefig(savename, bbox_inches='tight', dpi=100, decorate=True)
    plt.close()


def make_plots(df, variable, obs_variable, startdate, enddate, region,vmin,vmax,out_name):
    if startdate == None and enddate == None:
        for t in df.time.unique():
            date = pd.Timestamp(t)
            print(
                ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
            print('Creating Plot:', obs_variable, 'at time:', date)
            print(
                ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
            odf = df.loc[df.time ==
                          date, ['time', 'latitude', 'longitude', obs_variable, variable]]
            if ~odf.empty:
                make_spatial_bias_plot(
                    odf,
                    out_name,
                    vmin,
                    vmax,
                    col1=obs_variable,
                    col2=variable,
                    date=t,
                    region=region,
                    cmap='RdBu_r',
                    edgecolor='k',
                    linewidth=.8)
    else:
        sdate=pd.Timestamp(startdate)
        edate=pd.Timestamp(enddate)
        df_mean=df.groupby(['siteid'],as_index=False).mean()
        print(
                ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
        print('Creating Plot:', obs_variable, 'for period:', startdate, 'to ', enddate  )
        print(
                ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
        make_spatial_bias_plot(
                    df_mean,
                    out_name,
                    vmin,
                    vmax,
                    col1=obs_variable,
                    col2=variable,
                    date=edate,
                    region=region,
                    cmap='RdBu_r',
                    edgecolor='k',
                    linewidth=.8)

def get_df_region(obj, region):
    from monet.util.tools import get_giorgi_region_df as get_giorgi
    if region.lower() == 'domain':
        obj['GIORGI_ACRO'] = 'domain'
        return obj
    else:
        obj = get_giorgi(region)
        return obj.loc[obj.GIORGI_ACRO == region.upper()]


if __name__ == '__main__':

    parser = ArgumentParser(
        description='Make Spatial Plots for each time step or over period in files',
        formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '-p',
        '--paired_data',
        help='paired data input file names',
        type=str,
        nargs='+',
        required=True)
    parser.add_argument(
        '-s', '--species', nargs='+', help='Species', required=False, default=['aod_550nm'])
    parser.add_argument(
        '-b',
        '--subset_giorgi',
        help='Giorgi Region Subset true/false',
        type=bool,
        required=False,
        default=False)
    parser.add_argument(
        '-g',
        '--giorgi_region',
        help='GIORGI Region ACRONYMs NAU,SAU,AMZ,SSA,CAM,WNA,CNA,ENA,ALA,GRL,MED,NEU,WAF,EAF,SAF,SAH,SEA,EAS,SAS,CAS,TIB,NAS',
        required=False,
        default='domain')
    parser.add_argument(
        '-n',
        '--output_name',
        help='Spatial bias plot Output base name',
        type=str,
        required=False,
        default='AERONET_FV3CHEM')
    parser.add_argument(
        '-r',
        '--regulatory',
        help='boolean set to True fore 8-hrmax  or 24-ave NAAQS regulatory calcs',
        type=bool,
        required=False,
        default=False)
    parser.add_argument(
        '-sd',
        '--startdate',
        help='Startdate for bias plot statistics over a period YYYY-MM-DD HH:MM:SS',
        type=str,
        required=False,
        default=None)
    parser.add_argument(
        '-ed',
        '--enddate',
        help='Enddate for bias plot statisics over a period YYYY-MM-DD HH:MM:SS',
        type=str,
        required=False,
        default=None)
    parser.add_argument(
        '-miny', '--miny_scale', help='Set static min y-scale', type=float, required=False, default=None)
    parser.add_argument(
        '-maxy', '--maxy_scale', help='Set static max y-scale', type=float, required=False, default=None)
    args = parser.parse_args()

    paired_data = args.paired_data
    species     = args.species
    out_name    = args.output_name
    subset      = args.subset_giorgi
    region      = args.giorgi_region
    startdate   = args.startdate
    enddate     = args.enddate
    reg         = args.regulatory
    vmin        = args.miny_scale
    vmax        = args.maxy_scale


#load the paired dataframe 
    df = load_paired_data(paired_data)
    mapping_table = {'aod_550nm':'pm25aod550'}
    sub_map = {i: mapping_table[i] for i in species if i in mapping_table}
    if region is "domain":
     subset = False 
# subset  only the correct region
    #if subset is True:
    # df.query('giorgi_region == '+'"'+region+'"',inplace=True)

    #Loop through species
    for jj in species:
     df[jj] = np.where(df[jj]<=0, np.nan, df[jj]) #Replace all values < 0 with NaN
     df_drop=df.dropna(subset=[jj,sub_map.get(jj)]) #Drops all corresponding rows with obs species = NaN

#Converts OZONE, PM10, or PM2.5 dataframe to NAAQS regulatory values
     if jj == 'OZONE' and reg is True:
      df2 = make_8hr_regulatory(df_drop,[jj,sub_map.get(jj)]).rename(index=str,columns={jj+'_y':jj,sub_map.get(jj)+'_y':sub_map.get(jj)})
     elif jj == 'aod_550nm' and reg is True:
      df2 = make_24hr_regulatory(df_drop,[jj,sub_map.get(jj)]).rename(index=str,columns={jj+'_y':jj,sub_map.get(jj)+'_y':sub_map.get(jj)})
     elif jj == 'pm10_ugm3' and reg is True:
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
     #df2.rename(index=str,columns={"Q2": "RH_mod"},inplace=True)
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
       outname = "{}.{}.{}.{}.{}.{}".format(out_name,region, jj,startdatename, enddatename,'reg')
      if jj == 'PM2.5':
       outname = outname.replace('PM2.5','PM2P5')
      if region == 'domain':
       outname = outname.replace('domain','5X')
     else:
      dfnew = df2
      outname = "{}.{}.{}".format(out_name, region, jj)
      if reg is True:
       outname = "{}.{}.{}.{}".format(out_name,region, jj, 'reg')
      if jj == 'PM2.5':
       outname = outname.replace('PM2.5','PM2P5')
      if region == 'domain':
       outname = outname.replace('domain', '5X')

     dfnew_drop=dfnew.dropna(subset=[jj,sub_map.get(jj)])

     initial_datetime = dfnew_drop.time.min()
     # make the plots
     make_plots(dfnew_drop, sub_map.get(jj), jj, startdate, enddate, region,vmin,vmax,outname)
