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

import dask
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import cartopy.crs as ccrs
mpl.use('agg')

import pandas as pd
import seaborn as sns
import numpy as np
import monet
from monet.util.tools import calc_8hr_rolling_max,calc_24hr_ave,get_relhum

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


def map_projection(f):
    import cartopy.crs as ccrs
    if f.MAP_PROJ == 1:
        proj = ccrs.LambertConformal(
            central_longitude=f.CEN_LON, central_latitude=f.CEN_LAT)
    elif f.MAP_PROJ == 6:
        #Plate Carree is the equirectangular or equidistant cylindrical
        proj = ccrs.PlateCarree(
            central_longitude=f.CEN_LON)
    else:
        raise NotImplementedError('projection not supported')
    return proj


def chdir(fname):
    dir_path = os.path.dirname(os.path.realpath(fname))
    os.chdir(dir_path)
    return os.path.basename(fname)


def open_cmaq(finput):
    #from monetio.models import cmaq
    #f = cmaq.open_mfdataset(finput)
    import rapchem
    f = rapchem.open_mfdataset(finput)
    return f


def load_paired_data(fname):
    return pd.read_hdf(fname)

def make_spatial_plot(da, df, outname, proj, startdate, enddate,vmin,vmax,region='domain'): 
    cbar_kwargs = dict(aspect=30,shrink=.8)#dict(aspect=30)                       

    if region == 'domain':
     latmin= 20.0
     lonmin=-140.0
     latmax= 60.0
     lonmax=-60.0
    else:
     from monet.util.tools import get_epa_region_bounds as get_epa_bounds
     latmin,lonmin,latmax,lonmax,acro = get_epa_bounds(index=None,acronym=region)
    
    extent = [lonmin,lonmax,latmin,latmax]
    #change map_kwarg to map_kws
    ax = da.monet.quick_map(cbar_kwargs=cbar_kwargs, figsize=(15, 8), map_kws={'states': True, 'crs': proj,'extent':extent},robust=True,vmin=vmin,vmax=vmax,cmap=plt.cm.get_cmap('Spectral_r')) 
    plt.gcf().canvas.draw() 
    plt.tight_layout(pad=0)
    
    if startdate == None and enddate == None:
     date = pd.Timestamp(da.time.values) 
     dt = date - initial_datetime
     dtstr = str(dt.days * 24 + dt.seconds // 3600).zfill(3)
     plt.title(date.strftime('time=%Y/%m/%d %H:00 | RAP - AIRNOW '))
    else:
     plt.title('average time period | RAP - AIRNOW ')
     
    cbar = ax.figure.get_axes()[1] 
    if vmin == None and vmax == None:
     vmin, vmax = cbar.get_ybound()
     
    vars = df.keys() 
    varname = [x for x in vars if x not in ['latitude','longitude']][0] 
    ax.axes.scatter(df.longitude.values,df.latitude.values,s=25,c=df[varname],transform=ccrs.PlateCarree(),edgecolor='b',linewidth=.50,vmin=vmin,vmax=vmax,cmap=plt.cm.get_cmap('Spectral_r'))
    ax.axes.set_extent(extent,crs=ccrs.PlateCarree())
    
    if startdate == None and enddate == None:
     savename = "{}.{}.{}.jpg".format(outname,
                                     initial_datetime.strftime('sp'),
                                     dtstr)
    else:
     savename = "{}.{}.jpg".format(outname,
                                      'sp')
    print(savename)
    monet.plots.savefig(savename, bbox_inches='tight', dpi=100, decorate=True)
    plt.close()

def make_plots(finput, paired_data, variable, obs_variable, verbose, startdate, enddate, region,vmin,vmax, outname):
  if startdate == None and enddate == None:
    # open the files
    f = open_cmaq(finput)
    # get map projection
    proj = map_projection(f)
    if paired_data is not None:
        df = paired_data
    # loop over varaible list
    plots = []
    obj = f[variable]
        # loop over time
    for t in obj.time:
            date = pd.Timestamp(t.values)
            print(
                ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
            print('Creating Plot:', obs_variable, 'at time:', date)
            print(
                ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
            odf = df.loc[df.time == pd.Timestamp(t.values),['latitude','longitude',obs_variable]]
            make_spatial_plot(obj.sel(time=t), odf, outname, proj, startdate, enddate,vmin,vmax, region=region)
  else:
    # open the files
    f = open_cmaq(finput)
    # get map projection
    proj = map_projection(f)
    if paired_data is not None:
        df = paired_data
    # loop over varaible list
    plots = []
#    for index, var in enumerate(variable):
    obj = f[variable]
    sdate=pd.Timestamp(startdate)
    edate=pd.Timestamp(enddate)
    df_mean=df.groupby(['siteid'],as_index=False).mean()
    print(
        ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
    print('Creating Plot:', obs_variable, 'for period:', startdate, 'to ', enddate  )
    print(
        ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
    odf_mean = df_mean[['latitude','longitude',obs_variable]]
    mod_slice = obj.sel(time=slice(startdate,enddate)) 
    mod_mean = mod_slice.mean(dim='time')
    make_spatial_plot(mod_mean, odf_mean, outname, proj, startdate, enddate,vmin,vmax, region=region)

if __name__ == '__main__':

    parser = ArgumentParser(description='Make Spatial Plots for each time step in files',
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '-f', '--files', help='input model file names', nargs='+', type=str, required=True)
    parser.add_argument(
        '-p', '--paired_data', help='associated paired data input file name', type=str, required=True)
    parser.add_argument(
        '-s', '--species', nargs='+', help='species to plot', type=str,required=False, default=['OZONE'])
    parser.add_argument(
        '-v', '--verbose', help='print debugging information', action='store_true', required=False)
    parser.add_argument(
        '-b', '--subset_epa', help='EPA Region Subset true/false',type=bool, required=False, default=False)
    parser.add_argument(
        '-e', '--epa_region', help='EPA Region Acronyms', type=str, required=False, default='domain')
    parser.add_argument(
        '-n', '--output_name', help='Output base name',type=str, required=False, default='RAP_AIRNOW')
    parser.add_argument(
        '-sup', '--suppress_xwindow', help='Suppress X Window', action='store_true', required=False)
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

    finput      = args.files
    verbose     = args.verbose
    paired_data = args.paired_data
    species     = args.species
    out_name    = args.output_name
    subset      = args.subset_epa
    region      = args.epa_region
    startdate   = args.startdate
    enddate     = args.enddate
    reg         = args.regulatory
    vmin        = args.miny_scale
    vmax        = args.maxy_scale

    #load the paired dataframe
    
    df = load_paired_data(paired_data)
    #changed PM25_TOT to PM25
    mapping_table =  {
        'OZONE': 'o3',
        'PM2.5': 'PM25',
        'CO': 'CO',
        'NOX': 'NOx',
        'SO2': 'SO2',
        'NO': 'NO',
        'NO2': 'NO2',
    }
    sub_map = {i: mapping_table[i] for i in species if i in mapping_table}
    
    if region is "domain":
     subset = False 
# subset  only the correct region
    if subset is True:
     df.query('epa_region == '+'"'+region+'"',inplace=True)

    #Loop through species
    for jj in species:
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
      df2.loc[:,'WS']=df2.loc[:,'WS']*0.514  #convert obs knots-->m/s ***Conform to model units for overlay ***
      df2.query('WS > 0.2',inplace=True)  #Filter out calm WS obs (< 0.2 m/s), should not be trusted--creates artificially larger postive  model bias
     elif jj == 'BARPR':
      #df2.loc[:,'PRSFC']=df2.loc[:,'PRSFC']*0.01 #convert model Pascals-->millibars
      df2.loc[:,'BARPR']=df2.loc[:,'BARPR']/0.01 #convert obs millibars-->Pascals ***Conform to model units for overlay ***
     elif jj == 'PRECIP':
      df2.loc[:,'PRECIP']=df2.loc[:,'PRECIP']*0.1 #convert obs mm-->cm
     elif jj == 'TEMP':
      #df2.loc[:,'TEMP2'] = df2.loc[:,'TEMP2']-273.16 #convert model K-->C
      df2.loc[:,'TEMP'] = df2.loc[:,'TEMP']+273.16 #convert obs C-->K ***Conform to model units for overlay ***
     elif jj == 'RHUM':
     #convert model mixing ratio to relative humidity
      df2.loc[:,'Q2'] = get_relhum(df2.loc[:,'TEMP2'],df2.loc[:,'PRSFC'],df2.loc[:,'Q2'])  # *** Currently not supported for spatial overlay ***
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
      outname = "{}.{}.{}".format(out_name,region, jj)
      if reg is True:
       outname = "{}.{}.{}.{}".format(out_name,region, jj, 'reg')
      if jj == 'PM2.5':
       outname = outname.replace('PM2.5','PM2P5')
      if region == 'domain':
       outname = outname.replace('domain','5X')

     dfnew_drop=dfnew.dropna(subset=[jj,sub_map.get(jj)])

     initial_datetime = dfnew_drop.time.min()
     # make the plots



     make_plots(finput, dfnew_drop, sub_map.get(jj),
               jj, verbose, startdate, enddate,region,vmin,vmax, outname)
