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
Simple utility to make spatial plots from the NAQFC forecast
'''

initial_datetime = None

def map_projection(f):
    import cartopy.crs as ccrs
    proj = ccrs.LambertConformal(
        central_longitude=f.XCENT, central_latitude=f.YCENT)
    return proj


def chdir(fname):
    dir_path = os.path.dirname(os.path.realpath(fname))
    os.chdir(dir_path)
    return os.path.basename(fname)


def open_cmaq(finput):
    from monet.models import cmaq
    f = cmaq.open_mfdataset(finput)
    return f

def make_spatial_plot(da, outname, proj, startdate, enddate, count,vmin,vmax,region='domain'): 
    cbar_kwargs = dict(aspect=30,shrink=.8)#dict(aspect=30)                       

    if region == 'domain':
     latmin= 25.0
     lonmin=-130.0
     latmax= 48.0
     lonmax=-70.0
    else:
     from monet.util.tools import get_epa_region_bounds as get_epa_bounds
     latmin,lonmin,latmax,lonmax,acro = get_epa_bounds(index=None,acronym=region)
    
    extent = [lonmin,lonmax,latmin,latmax]
    ax = da.monet.quick_map(cbar_kwargs=cbar_kwargs, figsize=(15, 8), map_kwarg={'states': True, 'crs': proj,'extent':extent},robust=True,vmin=vmin,vmax=vmax,cmap=plt.cm.get_cmap('Spectral_r')) 
    plt.gcf().canvas.draw() 
    plt.tight_layout(pad=0)
    if startdate == None and enddate == None:
     date = pd.Timestamp(da.time.values )
     initial_datetime=date
#     dt = date - initial_datetime 
#     dtstr = str(dt.days * 24 + dt.seconds // 3600).zfill(3)
     dt = count
     dtstr = str(dt).zfill(3)
     plt.title(date.strftime('time=%Y/%m/%d %H:00 | CMAQ '))
    else:
     plt.title('average time period | CMAQ')

    cbar = ax.figure.get_axes()[1] 
    if vmin == None and vmax == None:
     vmin, vmax = cbar.get_ybound() 
#    vars = df.keys() 
#    varname = [x for x in vars if x not in ['latitude','longitude']][0] 
#    ax.scatter(df.longitude.values,df.latitude.values,s=25,c=df[varname],transform=ccrs.PlateCarree(),edgecolor='b',linewidth=.50,vmin=vmin,vmax=vmax,cmap=plt.cm.get_cmap('Spectral_r'))
    ax.set_extent(extent,crs=ccrs.PlateCarree())
    
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

def make_plots(finput, variable, verbose, startdate, enddate, region, vmin,vmax, outname):
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
    count = 0
    for t in obj.time:
            date = pd.Timestamp(t.values)
            print(
                ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
            print('Creating Plot:', variable, 'at time:', date)
            print(
                ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
            make_spatial_plot(obj.sel(time=t), outname, proj, startdate, enddate, count,vmin,vmax,region=region)
            count=count+1
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
    print(
        ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
    print('Creating Plot:', variable, 'for period:', startdate, 'to ', enddate  )
    print(
        ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
    mod_slice = obj.sel(time=slice(startdate,enddate)) 
    mod_mean = mod_slice.mean(dim='time')
    make_spatial_plot(mod_mean, outname, proj, startdate, enddate,count,vmin,vmax, region=region)

if __name__ == '__main__':

    parser = ArgumentParser(description='Make Spatial Plots for each time step in files',
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '-f', '--files', help='input model file names', nargs='+', type=str, required=True)
    parser.add_argument(
        '-p', '--paired_data', help='associated paired data input file name', type=str, required=False, default=None)
    parser.add_argument(
        '-s', '--species', nargs='+', help='model species to plot', type=str,required=False, default=['O3'])
    parser.add_argument(
        '-v', '--verbose', help='print debugging information', action='store_true', required=False)
    parser.add_argument(
        '-b', '--subset_epa', help='EPA Region Subset true/false',type=bool, required=False, default=False)
    parser.add_argument(
        '-e', '--epa_region', help='EPA Region Acronyms', type=str, required=False, default='domain')
    parser.add_argument(
        '-n', '--output_name', help='Output base name',type=str, required=False, default='CMAQ_AIRNOW')
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


    if region is "domain":
     subset = False 
    #Loop through species
    for jj in species:
    #subset for period, or use output frequency
     if startdate != None and enddate != None:
      import datetime
      startdatename_obj = datetime.datetime.strptime(startdate, '%Y-%m-%d %H:%M:%S')
      enddatename_obj   = datetime.datetime.strptime(enddate, '%Y-%m-%d %H:%M:%S')
      startdatename = str(datetime.datetime.strftime(startdatename_obj,'%Y-%m-%d_%H'))
      enddatename = str(datetime.datetime.strftime(enddatename_obj,'%Y-%m-%d_%H'))
      outname = "{}.{}.{}.{}.{}".format(out_name, region, jj, startdatename, enddatename)
      if jj == 'PM25_TOT':
       outname = outname.replace('PM25_TOT','PM2P5')
      if region == 'domain':
       outname = outname.replace('domain','5X')
     else:
      outname = "{}.{}.{}".format(out_name,region, jj)
      if jj == 'PM2.5':
       outname = outname.replace('PM25_TOT','PM2P5')
      if region == 'domain':
       outname = outname.replace('domain','5X')

     make_plots(finput, jj, verbose, startdate, enddate,region,vmin,vmax, outname)
