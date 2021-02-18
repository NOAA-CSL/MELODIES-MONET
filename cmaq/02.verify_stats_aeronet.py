#!/usr/bin/env python

__author__  = 'Patrick Campbell'
__email__   = 'patrick.c.campbell@noaa.gov'
__license__ = 'GPL'



#Simple MONET utility to calculate statistics from paired hdf file

import os
from glob import glob
import sys

import subprocess
from distutils.spawn import find_executable
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

import monet  
from monet.util.tools import calc_8hr_rolling_max,calc_24hr_ave,get_relhum
from monet.util.mystats import NO,NP,NOP,MO,MP,MdnO,MdnP,STDO,STDP,MB,WDMB_m,NMB,WDNMB_m,NMB_ABS,NME_m,NME_m_ABS,WDME_m,RMSE,WDRMSE_m,IOA_m,WDIOA_m,R2 
import pandas as pd
import numpy as np
from numpy import sqrt
import dask.dataframe as dd

#Define all statistics and statistical plots desired

def  calc_r2(obs,mod):
     """ Coefficient of Determination (unit squared) """
     return R2(obs,mod,axis=None)

def  calc_ioa_wd(obs,mod):
     """ Wind Direction Index of Agreement """
     return WDIOA_m(obs,mod,axis=0)

def  calc_ioa(obs,mod):
     """ Index of Agreement """
     return IOA_m(obs,mod,axis=0)

def  calc_rmse_wd(obs,mod):
     """ Wind Direction Root  Mean Square Error """
     return WDRMSE_m(obs,mod,axis=0)

def  calc_rmse(obs,mod):
     """ Root  Mean Square Error """
     return RMSE(obs,mod,axis=0)

def  calc_wdme(obs,mod):
     """ Wind Direction Mean Gross Error """
     return WDME_m(obs,mod,axis=0)

def  calc_nme_temp(obs,mod):
     """ Temperature (C) Normalized Mean Error (%)"""
     return NME_m_ABS(obs,mod,axis=0)

def  calc_nme(obs,mod):
     """ Normalized Mean Error (%)"""
     return NME_m(obs,mod,axis=0)

def  calc_nmb_wd(obs,mod):
     """ Wind Direction Temperature (C) Normalized Mean Bias (%)"""
     return WDNMB_m(obs,mod,axis=0)

def  calc_nmb_temp(obs,mod):
     """ Temperature (C) Normalized Mean Bias (%)"""
     return NMB_ABS(obs,mod,axis=0)

def  calc_nmb(obs,mod):
     """ Normalized Mean Bias (%)"""
     return NMB(obs,mod,axis=0)

def  calc_mb_wd(obs,mod):
     """ Wind Direction Mean Bias """
     return WDMB_m(obs,mod,axis=0)

def  calc_mb(obs,mod):
     """ Mean Bias """
     return MB(obs,mod,axis=0)

def  calc_stdp(obs,mod):
     """ Standard deviation of Predictions """
     return STDP(obs,mod,axis=0)

def  calc_stdo(obs,mod):
     """ Standard deviation of Observations """
     return STDO(obs,mod,axis=0)

def  calc_MdnP(obs,mod):
     """ Median Predictions (obs unit) """
     return MdnP(obs,mod,axis=0)

def  calc_MdnO(obs,mod):
     """ Median Observations (model unit) """
     return MdnO(obs,mod,axis=0)

def  calc_MP(obs,mod):
     """ Mean Predictions (model unit) """
     return MP(obs,mod,axis=0)

def  calc_MO(obs,mod):
     """ Mean Observations (obs unit) """
     return MO(obs,mod,axis=0)

def  calc_NOP(obs,mod):
     """ N Observations/Prediction Pairs (#)"""
     return NOP(obs,mod,axis=0)

def  calc_NP(obs,mod):
     """ N Predictions (#) """
     return NP(obs,mod,axis=0)

def  calc_NO(obs,mod):
     """ N Observations (#) """
     return NO(obs,mod,axis=0)

def  make_24hr_regulatory(df,col=None):
     """ Make 24-hour averages """
     return calc_24hr_ave(df,col)

def  make_8hr_regulatory(df,col=None):
     """ Make 8-hour rolling average daily """
     return calc_8hr_rolling_max(df,col,window=8)


if __name__ == '__main__':

    parser = ArgumentParser(description='calculates statistics from paired file', formatter_class=ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('-p',   '--paired_data',       help='string input paired file directory/names', type=str,nargs='+', required=True)
    parser.add_argument('-r',   '--regulatory',  help='boolean set to True fore 8-hrmax  or 24-ave NAAQS regulatory calcs', type=bool, required=False, default=False)
    parser.add_argument('-sd',  '--startdate',   help='string start date to isolate periods for statistics YYYY-MM-DD HH:MM:SS', type=str, required=False, default=None)
    parser.add_argument('-ed',  '--enddate',     help='string end date to isolate periods for statistics YYYY-MM-DD HH:MM:SS', type=str, required=False, default=None)
    parser.add_argument('-s',   '--species',     help='string/list input for obs species-variables to create stats',type=str,nargs='+', required=False, default=['aod_550nm'])
    parser.add_argument('-b',   '--subset_giorgi',  help='boolean set to True for subsetting by GIORGI region', type=bool, required=False, default=False)
    parser.add_argument('-g',   '--giorgi_regions', help='string/list input for set GIORGI regions',type=str,nargs='+', required=False, default=['R1'])
    parser.add_argument('-v',   '--verbose',     help='print debugging information', action='store_true', required=False)
    args = parser.parse_args()

    finput       = args.paired_data
    reg          = args.regulatory
    startdate    = args.startdate
    enddate      = args.enddate
    species      = args.species
    subset_giorgi   = args.subset_giorgi
    giorgi_regions  = args.giorgi_regions
    verbose      = args.verbose

    for ee in giorgi_regions:
       df = dd.read_hdf(finput, '/*').compute()
       mapping_table = {'aod_550nm':'pm25aod550'}
       sub_map = {i: mapping_table[i] for i in species if i in mapping_table}   
#subsetting data for dates, regulatory calc, and/or giorgi regions    
       if startdate != None and enddate != None:
          mask = (df['time'] >= startdate) & (df['time'] <= enddate)
          df =df.loc[mask]
          import datetime
          startdatename_obj = datetime.datetime.strptime(startdate, '%Y-%m-%d %H:%M:%S')
          enddatename_obj   = datetime.datetime.strptime(enddate, '%Y-%m-%d %H:%M:%S')
          startdatename = str(datetime.datetime.strftime(startdatename_obj,'%Y-%m-%d_%H'))
          enddatename = str(datetime.datetime.strftime(enddatename_obj,'%Y-%m-%d_%H'))
       else:
          startdatename='Entire'
          enddatename  ='Period'

       if subset_giorgi is True:
          #df.query('giorgi_region == '+'"'+ee+'"',inplace=True)
          from monet.util.tools import get_giorgi_region_bounds as get_giorgi_bounds
          latmin,lonmin,latmax,lonmax,acro = get_giorgi_bounds(index=None,acronym=ee)
          df = df[(df['latitude'] >= latmin) & (df['latitude'] <= latmax)]
          df = df[(df['longitude'] >= lonmin) & (df['longitude'] <= lonmax)]  

       if reg is True and subset_giorgi is False:
          stats=open(finput[0].replace('.hdf','_')+startdatename+'_'+enddatename+'_reg_stats_domain.txt','w')
       elif reg is True and subset_giorgi is True:
          stats=open(finput[0].replace('.hdf','_')+startdatename+'_'+enddatename+'_reg_stats_'+ee+'.txt','w')
       elif reg is False and subset_giorgi is True:
          stats=open(finput[0].replace('.hdf','_')+startdatename+'_'+enddatename+'_stats_'+ee+'.txt','w')
       else:
          stats=open(finput[0].replace('.hdf','_')+startdatename+'_'+enddatename+'_stats_domain.txt','w')          
       
#Converts OZONE, PM10, or PM2.5 dataframe to NAAQS regulatory values
       for jj in species: 
        df[jj] = np.where(df[jj]<=0, np.nan, df[jj]) #Replace all values < 0 with NaN
       	df_drop=df.dropna(subset=[jj,sub_map.get(jj)]) #Drops all corresponding rows with obs species = NaN        
       
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
         df2.query('WS > 0.2',inplace=True)  #Filter out calm WS obs (< 0.2 m/s), calm obs winds should not be trusted--creates artificially larger postive  model bias
        elif jj == 'BARPR':
         df2.loc[:,'PRSFC']=df2.loc[:,'PRSFC']*0.01 #convert model Pascals-->millibars
        elif jj == 'PRECIP':
         df2.loc[:,'PRECIP']=df2.loc[:,'PRECIP']*0.1 #convert obs mm-->cm
        elif jj == 'TEMP':
         df2.loc[:,'TEMP2'] = df2.loc[:,'TEMP2']-273.16 #convert model K-->C
        elif jj == 'RHUM':
         #convert model mixing ratio to relative humidity
         df2.loc[:,'Q2'] = get_relhum(df2.loc[:,'TEMP2'],df2.loc[:,'PRSFC'],df2.loc[:,'Q2'])
         df2.rename(index=str,columns={"Q2": "RH_mod"},inplace=True)
        elif jj == 'CO':
         df2.loc[:,'CO']=df2.loc[:,'CO']*1000.0 #convert obs ppm-->ppb
        else:
         df2=df2  
#Calculates average statistics over entire file time
        if reg is True and subset_giorgi is False:
       	 stats=open(finput[0].replace('.hdf','_')+startdatename+'_'+enddatename+'_reg_stats_domain.txt','a')
        elif reg is True and subset_giorgi is True:
         stats=open(finput[0].replace('.hdf','_')+startdatename+'_'+enddatename+'_reg_stats_'+ee+'.txt','a')
        elif reg is False and subset_giorgi is True:
         stats=open(finput[0].replace('.hdf','_')+startdatename+'_'+enddatename+'_stats_'+ee+'.txt','a')
        else:
       	 stats=open(finput[0].replace('.hdf','_')+startdatename+'_'+enddatename+'_stats_domain.txt','a')
         
        print('---------------------------------')
        print('Statistics of ',jj,' and ',sub_map.get(jj),' pair over ', startdatename,' to ', enddatename)
        stats.write('Statistics of '+jj+' and '+sub_map.get(jj)+' pair over file period '+startdatename+ ' to '+enddatename+'\n')
        stats.write('---------------------------------'+'\n')

        obs_stats = df2[jj]
        if jj == 'RHUM':
         mod_stats = df2['RH_mod']
        else:
         mod_stats = df2[sub_map.get(jj)]
        no_stats=calc_NO(obs_stats,mod_stats)
        print('Number of',jj,' Observations = ', no_stats)
        stats.write('Number of '+jj+' Observations = '+str(no_stats)+'\n')
        
        np_stats=calc_NP(obs_stats,mod_stats)
        print('Number of',sub_map.get(jj),' Predictions = ', np_stats)
        stats.write('Number of '+sub_map.get(jj)+' Predictions = '+str(np_stats)+'\n')        

        nop_stats=calc_NOP(obs_stats,mod_stats)
        print('Number of',jj,'/',sub_map.get(jj),' Observations/Prediction Pairs (#) = ', nop_stats)
        stats.write('Number of '+jj+'/'+sub_map.get(jj)+' Observations/Prediction Pairs (#) = '+str(nop_stats)+'\n')        

        mo_stats=calc_MO(obs_stats,mod_stats)
        print('Mean of',jj,' Observations = ', "{:8.2f}".format(mo_stats))
        stats.write('Mean of '+jj+' Observations = '+str("{:8.2f}".format(mo_stats))+'\n') 
       
        mp_stats=calc_MP(obs_stats,mod_stats)
        print('Mean of',sub_map.get(jj),' Predictions = ', "{:8.2f}".format(mp_stats))
        stats.write('Mean of '+sub_map.get(jj)+' Predictions = '+str("{:8.2f}".format(mp_stats))+'\n')       
 
        mdno_stats=calc_MdnO(obs_stats,mod_stats)
        print('Median of',jj,' Observations = ', "{:8.2f}".format(mdno_stats))
        stats.write('Median of '+jj+' Observations = '+str("{:8.2f}".format(mdno_stats))+'\n')        

        mdnp_stats=calc_MdnP(obs_stats,mod_stats)
        print('Median of',sub_map.get(jj),' Predictions = ', "{:8.2f}".format(mdnp_stats))
        stats.write('Median of '+sub_map.get(jj)+' Predictions = '+str("{:8.2f}".format(mdnp_stats))+'\n')

        stdo_stats=calc_stdo(obs_stats,mod_stats)
        print('Standard deviation of',jj,' Observations = ', "{:8.2f}".format(stdo_stats))
        stats.write('Standard deviation of '+jj+' Observations = '+str("{:8.2f}".format(stdo_stats))+'\n')
       
        stdp_stats=calc_stdp(obs_stats,mod_stats)
        print('Standard deviation of',sub_map.get(jj),' Predictions = ', "{:8.2f}".format(stdp_stats))
        stats.write('Standard deviation of '+sub_map.get(jj)+' Predictions = '+str("{:8.2f}".format(stdp_stats))+'\n')       

        if jj == 'WD':
         mb_stats=calc_mb_wd(obs_stats,mod_stats)
         print('Wind Direction Mean Bias of ',sub_map.get(jj),'-',jj,' =  ', "{:8.2f}".format(mb_stats))
         stats.write('Wind Direction Mean Bias of '+sub_map.get(jj)+'-'+jj+' =  '+str("{:8.2f}".format(mb_stats))+'\n')
        else:
         mb_stats=calc_mb(obs_stats,mod_stats)
         print('Mean Bias of ',sub_map.get(jj),'-',jj,' =  ', "{:8.2f}".format(mb_stats))
         stats.write('Mean Bias of '+sub_map.get(jj)+'-'+jj+' =  '+str("{:8.2f}".format(mb_stats))+'\n')        
        
        if jj == 'TEMP':
         nmb_stats=calc_nmb_temp(obs_stats,mod_stats)
         print('Temperature Normalized Mean Bias (%) of ',sub_map.get(jj),'-',jj,' =  ', "{:8.2f}".format(nmb_stats))        
         stats.write('Temperature Normalized Mean Bias (%) of '+sub_map.get(jj)+'-'+jj+' =  '+str("{:8.2f}".format(nmb_stats))+'\n')
        elif jj == 'WD':
         nmb_stats=calc_nmb_wd(obs_stats,mod_stats)
         print('Wind Direction Normalized Mean Bias (%) of ',sub_map.get(jj),'-',jj,' =  ', "{:8.2f}".format(nmb_stats))
         stats.write('Wind Direction Normalized Mean Bias (%) of '+sub_map.get(jj)+'-'+jj+' =  '+str("{:8.2f}".format(nmb_stats))+'\n')
        else:
         nmb_stats=calc_nmb(obs_stats,mod_stats)
         print('Normalized Mean Bias (%) of ',sub_map.get(jj),'-',jj,' =  ', "{:8.2f}".format(nmb_stats))
         stats.write('Normalized Mean Bias (%) of '+sub_map.get(jj)+'-'+jj+' =  '+str("{:8.2f}".format(nmb_stats))+'\n')
        
        if jj == 'TEMP':
         nme_stats=calc_nme_temp(obs_stats,mod_stats)
         print('Temperature Normalized Mean Error (%) of ',sub_map.get(jj),'-',jj,' =  ', "{:8.2f}".format(nme_stats))
         stats.write('Temperature Normalized Mean Error (%) of '+sub_map.get(jj)+'-'+jj+' =  '+str("{:8.2f}".format(nme_stats))+'\n')
        elif jj == 'WD':
         wdme_stats=calc_wdme(obs_stats,mod_stats)
         print('Wind Direction Mean Gross Error of ',sub_map.get(jj),'-',jj,' =  ', "{:8.2f}".format(wdme_stats))
         stats.write('Wind Direction Mean Gross Error (%) of '+sub_map.get(jj)+'-'+jj+' =  '+str("{:8.2f}".format(wdme_stats))+'\n')
        else:
         nme_stats=calc_nme(obs_stats,mod_stats)
         print('Normalized Mean Error (%) of ',sub_map.get(jj),'-',jj,' =  ', "{:8.2f}".format(nme_stats))
         stats.write('Normalized Mean Error (%) of '+sub_map.get(jj)+'-'+jj+' =  '+str("{:8.2f}".format(nme_stats))+'\n')
        
        if jj == 'WD':
         rmse_stats=calc_rmse_wd(obs_stats,mod_stats)
         print('Wind Direction Root Mean Square Error of ',sub_map.get(jj),'-',jj,' =  ', "{:8.2f}".format(rmse_stats))
         stats.write('Wind Direction Root Mean Square Error of '+sub_map.get(jj)+'-'+jj+' =  '+str("{:8.2f}".format(rmse_stats))+'\n')
        else:
         rmse_stats=calc_rmse(obs_stats,mod_stats)
         print('Root Mean Square Error of ',sub_map.get(jj),'-',jj,' =  ', "{:8.2f}".format(rmse_stats))
         stats.write('Root Mean Square Error of '+sub_map.get(jj)+'-'+jj+' =  '+str("{:8.2f}".format(rmse_stats))+'\n')
        
        if jj == 'WD':
         ioa_stats=calc_ioa_wd(obs_stats,mod_stats)
         print('Wind Direction Index of Agreement of ',sub_map.get(jj),'-',jj,' =  ', "{:8.2f}".format(ioa_stats))
         stats.write('Wind Direction Index of Agreement of '+sub_map.get(jj)+'-'+jj+' =  '+str("{:8.2f}".format(ioa_stats))+'\n')
        else:
         ioa_stats=calc_ioa(obs_stats,mod_stats)
         print('Index of Agreement of ',sub_map.get(jj),'-',jj,' =  ', "{:8.2f}".format(ioa_stats))
         stats.write('Index of Agreement of '+sub_map.get(jj)+'-'+jj+' =  '+str("{:8.2f}".format(ioa_stats))+'\n')

        r_stats=sqrt(calc_r2(obs_stats,mod_stats))
        print('Pearsons Correlation Coefficient of ',sub_map.get(jj),'-',jj,' =  ', "{:8.2f}".format(r_stats))
        stats.write('Pearsons Correlation Coefficient of '+sub_map.get(jj)+'-'+jj+' =  '+str("{:8.2f}".format(r_stats))+'\n')

        print('Statistics done!')
        print('---------------------------------')
        stats.write('---------------------------------'+'\n')
        stats.close()

    sys.exit(0)
    


